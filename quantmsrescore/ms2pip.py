import itertools
import logging
import re
from collections import defaultdict
from itertools import chain

from duckdb.duckdb import identifier
from math import ceil
from pathlib import Path
from typing import Optional, Tuple, List, Union, Callable, Generator

import numpy as np
from ms2pip import correlate
from ms2pip._cython_modules import ms2pip_pyx
from ms2pip._utils.encoder import Encoder
from ms2pip._utils.ion_mobility import IonMobility
from ms2pip._utils.psm_input import read_psms
from ms2pip._utils.retention_time import RetentionTime
from ms2pip.constants import MODELS
from ms2pip.core import _Parallelized, _process_peptidoform
from ms2pip.exceptions import NoMatchingSpectraFound
from ms2pip.result import ProcessingResult, calculate_correlations
from ms2pip.spectrum import ObservedSpectrum
from ms2rescore.feature_generators import MS2PIPFeatureGenerator
from ms2rescore.feature_generators.base import FeatureGeneratorException
from ms2rescore.utils import infer_spectrum_path
from psm_utils import PSMList, PSM
from sdrf_pipelines.openms.openms import OpenMS

from quantmsrescore.constants import SUPPORTED_MODELS_MS2PIP
from quantmsrescore.exceptions import Ms2pipIncorrectModelException
import ms2pip.exceptions as exceptions

from quantmsrescore.openms import OpenMSHelper

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

class PatchParallelized(_Parallelized):
    """
    Extended version of _Parallelized that supports custom spectrum file reading using pyopenms instead of
    ms2rescore_rs.
    """

    def __init__(
            self,
            encoder,
            model=None,
            model_dir=None,
            ms2_tolerance=0.02,
            processes=None
    ):
        """
        Initialize with all original parameters plus a custom spectrum reader.

        Parameters
        ----------
        spectrum_reader : Callable
            Function that reads spectrum files and yields spectrum objects
        """
        super().__init__(
            encoder=encoder,
            model=model,
            model_dir=model_dir,
            ms2_tolerance=ms2_tolerance,
            processes=processes,
        )

    def process_spectra(
            self,
            psm_list,
            spectrum_file,
            spectrum_id_pattern,
            vector_file=False,
            annotations_only=False,
    ):
        """
        Override process_spectra to use our custom spectrum reader
        """
        # Validate runs and collections
        if not len(psm_list.collections) == 1 or not len(psm_list.runs) == 1:
            raise exceptions.InvalidInputError("PSMs should be for a single run and collection.")

        # Define our custom _process_spectra function that uses our reader
        # Use our custom function with the execute_in_pool method
        args = (
            spectrum_file,
            vector_file,
            self.encoder,
            self.model,
            self.ms2_tolerance,
            spectrum_id_pattern,
            annotations_only,
        )

        results = self._execute_in_pool(psm_list, _custom_process_spectra, args)

        # Validate number of results
        if not results:
            raise exceptions.NoMatchingSpectraFound(
                "No spectra matching spectrum IDs from PSM list could be found in provided file."
            )
        logging.debug(f"Gathered data for {len(results)} PSMs.")

        # Add XGBoost predictions if required
        if (
                not (vector_file or annotations_only)
                and "xgboost_model_files" in MODELS[self.model].keys()
        ):
            results = self._add_xgboost_predictions(results)

        return results

    def _execute_in_pool(self, psm_list: PSMList, func: Callable, args: tuple):
        """Execute function in multiprocessing pool."""

        def get_chunk_size(n_items, n_processes):
            """Get optimal chunk size for multiprocessing."""
            if n_items < 5000:
                return n_items
            else:
                max_chunk_size = 50000
                n_chunks = ceil(ceil(n_items / n_processes) / max_chunk_size) * n_processes
                return ceil(n_items / n_chunks)

        def to_chunks(_list, chunk_size):
            """Split _list into chunks of size chunk_size."""

            def _generate_chunks():
                for i in range(0, len(_list), chunk_size):
                    yield _list[i : i + chunk_size]

            _list = list(_list)
            return list(_generate_chunks())

        def _enumerated_psm_list_by_spectrum_id(psm_list, spectrum_ids_chunk):
            selected_indices = np.flatnonzero(np.isin(psm_list["spectrum_id"], spectrum_ids_chunk))
            return [(i, psm_list.psm_list[i]) for i in selected_indices]

        with self._get_pool() as pool:
            if not psm_list:
                logging.warning("No PSMs to process.")
                return []

            # Split PSMList into chunks
            if func == _custom_process_spectra:
                # Split by spectrum_id to keep PSMs for same spectrum together
                spectrum_ids = set(psm_list["spectrum_id"])
                chunk_size = get_chunk_size(len(spectrum_ids), pool._processes)
                chunks = [
                    _enumerated_psm_list_by_spectrum_id(psm_list, spectrum_ids_chunk)
                    for spectrum_ids_chunk in to_chunks(spectrum_ids, chunk_size)
                ]
            else:
                # Simple split by PSM
                chunk_size = get_chunk_size(len(psm_list), pool._processes)
                chunks = to_chunks(list(enumerate(psm_list)), chunk_size)

            logging.debug(f"Processing {len(chunks)} chunk(s) of ~{chunk_size} entries each.")

            # Add jobs to pool
            mp_results = []
            for psm_list_chunk in chunks:
                mp_results.append(pool.apply_async(func, args=(psm_list_chunk, *args)))

            # Gather results
            # results = [
            #     r.get()
            #     for r in track(
            #         mp_results,
            #         disable=len(chunks) == 1,
            #         description="Processing chunks...",
            #         transient=True,
            #         show_speed=False,
            #     )
            # ]
            results = [r.get() for r in mp_results]

        # Sort results by input order
        results = list(
            sorted(
                itertools.chain.from_iterable(results),
                key=lambda result: result.psm_index,
            )
        )

        return results

    def spectrum_reader(self, spec_file) -> List:
        logging.info("Reading with Pyopenms")


class MS2PIPAnnotator(MS2PIPFeatureGenerator):

    def __init__(
        self,
        *args,
        model: str = "HCD",
        ms2_tolerance: float = 0.02,
        spectrum_path: Optional[str] = None,
        spectrum_id_pattern: str = "(.*)",
        model_dir: Optional[str] = None,
        processes: int = 1,
        calibration_set_size: Optional[float] = 0.20,
        valid_correlations_size: Optional[float] = 0.80,
        correlation_threshold: Optional[float] = 0.6,
        higher_score_better: bool = True,
        annotated_ms_tolerance: Tuple[float, str] = (0.0, None),
        predicted_ms_tolerance: Tuple[float, str] = (0.0, None),
        **kwargs,
    ):
        super().__init__(
            args,
            model=model,
            ms2_tolerance=ms2_tolerance,
            spectrum_path=spectrum_path,
            spectrum_id_pattern=spectrum_id_pattern,
            model_dir=model_dir,
            processes=processes,
            kwargs=kwargs,
        )
        self._reported_tolerance: Tuple[float, str] = annotated_ms_tolerance
        self._predicted_tolerance: Tuple[float, str] = predicted_ms_tolerance
        self._calibration_set_size: float = calibration_set_size
        self._valid_correlations_size: float = valid_correlations_size
        self._correlation_threshold: float = correlation_threshold
        self._higher_score_better: bool = higher_score_better

    def add_features(self, psm_list: PSMList) -> None:
        """
        Add MS²PIP-derived features to PSMs.

        Parameters
        ----------
        psm_list
        PSMs to add features to.
        """
        logging.info("Adding MS²PIP-derived features to PSMs.")
        psm_dict = psm_list.get_psm_dict()
        current_run = 1
        total_runs = sum(len(runs) for runs in psm_dict.values())

        for runs in psm_dict.values():
            for run, psms in runs.items():
                logging.info(
                    f"Running MS²PIP {self.model} for PSMs from run ({current_run}/{total_runs}) `{run}`..."
                )
                psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
                spectrum_filename = infer_spectrum_path(self.spectrum_path, run)
                logging.debug(f"Using spectrum file `{spectrum_filename}`")
                try:
                    ms2pip_results = self.custom_correlate(
                        psms=psm_list_run,
                        spectrum_file=str(spectrum_filename),
                        spectrum_id_pattern=self.spectrum_id_pattern,
                        model=self.model,
                        ms2_tolerance=self.ms2_tolerance,
                        compute_correlations=True,
                        model_dir=self.model_dir,
                        processes=self.processes,
                    )
                except NoMatchingSpectraFound as e:
                    raise FeatureGeneratorException(
                        f"Could not find any matching spectra for PSMs from run `{run}`. "
                        "Please check that the `spectrum_id_pattern` and `psm_id_pattern` "
                        "options are configured correctly. See "
                        "https://ms2rescore.readthedocs.io/en/latest/userguide/configuration/#mapping-psms-to-spectra"
                        " for more information."
                    ) from e
                valid_correlation = self._validate_scores(
                    ms2pip_results=ms2pip_results,
                    calibration_set_size=self._calibration_set_size,
                    valid_correlations_size=self._valid_correlations_size,
                    correlation_threshold=self._correlation_threshold,
                    higher_score_better=self._higher_score_better,
                )

                if not valid_correlation:
                    logging.error(
                        "The number of valid correlations doesn't exceed the threshold for current the calibration set."
                        "Please try a different model or adjust the valid_correlations_size or calibration_set_size."
                    )
                    raise Ms2pipIncorrectModelException(
                        message="The number of valid correlations doesn't exceed the threshold for current the "
                                "calibration set. Please try a different model or adjust the valid_correlations_size "
                                "or calibration_set_size.",
                        model=self.model,
                    )
                self._calculate_features(psm_list_run, ms2pip_results)
                current_run += 1

    @staticmethod
    def _validate_scores(
        ms2pip_results, calibration_set_size, valid_correlations_size, correlation_threshold, higher_score_better
    ) -> bool:
        """
        Validate MS²PIP results based on score and correlation criteria.

        This method checks if the MS²PIP results meet the specified correlation
        threshold and score criteria. It first filters out decoy PSMs, sorts the
        results based on the PSM score, and selects a calibration set. The method
        then verifies if at least 80% of the calibration set has a correlation
        above the given threshold.

        Parameters
        ----------
        ms2pip_results : list
            List of MS²PIP results to validate.
        calibration_set_size : float
            Fraction of the results to use for calibration.
        valid_correlations_size: float
            Fraction of the valid PSM.
        correlation_threshold : float
            Minimum correlation value required for a result to be considered valid.
        higher_score_better : bool
            Indicates if a higher PSM score is considered better.

        Returns
        -------
        bool
            True if the results are valid based on the criteria, False otherwise.
        """
        if not ms2pip_results:
            return False

        ms2pip_results_copy = (
            ms2pip_results.copy()
        )  # Copy ms2pip results to avoid modifying the original list

        # Select only PSMs that are target and not decoys
        ms2pip_results_copy = [result for result in ms2pip_results_copy if not result.psm.is_decoy]
        # Sort ms2pip results by PSM score and lower score is better
        ms2pip_results_copy.sort(key=lambda x: x.psm.score, reverse=higher_score_better)

        # Get a calibration set, the % of psms to be used for calibrarion is defined by calibration_set_size
        calibration_set = ms2pip_results_copy[
            : int(len(ms2pip_results_copy) * calibration_set_size)
        ]

        # Select the results with correlation above the threshold
        valid_correlation = [
            psm for psm in calibration_set if psm.correlation >= correlation_threshold
        ]

        logging.info(
            f"The number of valid correlations is {int(len(valid_correlation)/len(calibration_set)*100)}% of the "
            f"calibration set top {calibration_set_size*100}% PSMs"
        )

        # If the number of valid correlations is less than 80% of the calibration set, return False
        if len(valid_correlation) < len(calibration_set) * valid_correlations_size:
            return False

        return True

    def _find_best_ms2pip_model(
        self, batch_psms: PSMList, knwon_fragmentation: Optional[str] = None
    ) -> Tuple[str, float, float]:
        """
        Find the best MS²PIP model for a batch of PSMs.

        This method finds the best MS²PIP model for a batch of PSMs by
        comparing the correlation of the PSMs with the different models.

        Parameters
        ----------
        batch_psms : list
            List of PSMs to find the best model for.

        Returns
        -------
        Tuple
            Tuple containing the best model and the correlation value.
        """
        best_model = None
        best_correlation = 0

        filtered_models = SUPPORTED_MODELS_MS2PIP

        if knwon_fragmentation:
            filtered_models = {
                knwon_fragmentation: SUPPORTED_MODELS_MS2PIP.get(knwon_fragmentation)
            }

        self.ms2_tolerance = self.choose_best_ms2pip_tolerance(
            ms2_tolerance=self.ms2_tolerance, reported_tolerance=self._reported_tolerance
            , predicted_tolerance=self._predicted_tolerance
        )

        for fragment_types in filtered_models:
            for model in filtered_models[fragment_types]:
                logging.info(f"Running MS²PIP for model `{model}`...")
                ms2pip_results = correlate(
                    psms=batch_psms,
                    spectrum_file=self.spectrum_path,
                    spectrum_id_pattern=self.spectrum_id_pattern,
                    model=model,
                    ms2_tolerance=self.ms2_tolerance,
                    compute_correlations=True,
                    model_dir=self.model_dir,
                    processes=self.processes,
                )
                correlation = self._calculate_correlation(ms2pip_results)
                if correlation > best_correlation and correlation >= 0.4:
                    best_model = model
                    best_correlation = correlation

        return best_model, best_correlation, self.ms2_tolerance

    @staticmethod
    def _calculate_correlation(ms2pip_results: List[ProcessingResult]) -> float:
        """
        Calculate the average correlation from MS²PIP results.

        This method computes the average correlation score from a list of
        MS²PIP results, where each result contains a correlation attribute.

        Parameters
        ----------
        ms2pip_results : list
            List of MS²PIP results, each containing a correlation score.

        Returns
        -------
        float
            The average correlation score of the provided MS²PIP results.
        """
        total_correlation = sum([psm.correlation for psm in ms2pip_results])
        return total_correlation / len(ms2pip_results)

    def choose_best_ms2pip_tolerance(
        self,
        ms2_tolerance: float,
        reported_tolerance: Tuple[float, str],
        predicted_tolerance: Tuple[float, str] = None,
    ) -> float:
        """
        Determine the best MS²PIP tolerance to use by comparing tolerances in Da.

        Logic:
          - If reported tolerance is None: use the command line tolerance
          - If reported tolerance is in ppm: try to use the predicted Da equivalent if available
            and if it's not too different from the command line value
          - If command line tolerance is less restrictive (higher): use command line value
          - If command line tolerance is more restrictive (lower):
          - Use reported tolerance if command line is at least 10% of reported value
          - Otherwise use command line tolerance
          - Log all decisions for transparency

          Parameters
          ----------
          ms2_tolerance : float
              The MS²PIP tolerance specified in the command line (in Da).
          reported_tolerance : tuple(float, str)
              The tolerance reported in the idXML as (value, unit).
              Unit can be 'Da', 'ppm', or None.
          predicted_tolerance : tuple(float, str), optional
             The tolerance predicted in Da if the reported one is in ppm.
             Format is (value, unit) where unit should be 'Da' if provided.

          Returns
          -------
          float
             The best tolerance value to use (in Da).
        """

        # Case 1: No reported tolerance
        if reported_tolerance[1] is None:
            logging.info(
                f"No MS²PIP tolerance reported in the idXML. Using command line value ({ms2_tolerance} Da)."
            )
            return ms2_tolerance

        # Case 2: Reported tolerance is in ppm
        if reported_tolerance[1] == "ppm":
            if (
                predicted_tolerance is not None
                and predicted_tolerance[1] == "Da"
                and ms2_tolerance < predicted_tolerance[0]
                and (predicted_tolerance[0]/ms2_tolerance) > 0.1
            ):
                logging.warning(
                    f"Reported MS²PIP tolerance is in ppm. Using the predicted Da equivalent: "
                    f"{predicted_tolerance[0]} Da (instead of command line value: {ms2_tolerance} Da)."
                )
                return predicted_tolerance[0]
            else:
                logging.warning(
                    f"Reported MS²PIP tolerance is in ppm and no suitable Da prediction available. "
                    f"Using command line value ({ms2_tolerance} Da)."
                )
                return ms2_tolerance

        # Command line value is more restrictive (lower)
        if ms2_tolerance < reported_tolerance[0]:
            ratio = ms2_tolerance / reported_tolerance[0]
            if ratio > 0.1:
                logging.warning(
                    f"Command line MS²PIP tolerance ({ms2_tolerance} Da) is more restrictive "
                    f"by {(1 - ratio) * 100:.1f}% than reported value ({reported_tolerance[0]} Da). "
                    f"Using reported value for better model compatibility."
                )
                return reported_tolerance[0]

            else:
                logging.warning(
                    f"Command line MS²PIP tolerance ({ms2_tolerance} Da) is significantly more restrictive "
                    f"({ratio * 100:.1f}% of reported value: {reported_tolerance[0]} Da). "
                    f"Using command line value as specified."
                )
                return ms2_tolerance

        # Values are the same
        logging.info(
            f"MS²PIP tolerance in command line ({ms2_tolerance} Da) will continue be used instead of any reported idXML value."
        )
        return ms2_tolerance

    def custom_correlate(
            self,
            psms: Union[PSMList, str, Path],
            spectrum_file: Union[str, Path],
            psm_filetype: Optional[str] = None,
            spectrum_id_pattern: Optional[str] = None,
            compute_correlations: bool = False,
            add_retention_time: bool = False,
            add_ion_mobility: bool = False,
            model: Optional[str] = "HCD",
            model_dir: Optional[Union[str, Path]] = None,
            ms2_tolerance: float = 0.02,
            processes: Optional[int] = None,
    ) -> List[ProcessingResult]:
        """
        Custom implementation of correlate that uses our custom spectrum reader.
        """
        psm_list = read_psms(psms, filetype=psm_filetype)
        spectrum_id_pattern = spectrum_id_pattern if spectrum_id_pattern else "(.*)"

        if add_retention_time:
            logging.info("Adding retention time predictions")
            rt_predictor = RetentionTime(processes=processes)
            rt_predictor.add_rt_predictions(psm_list)

        if add_ion_mobility:
            logging.info("Adding ion mobility predictions")
            im_predictor = IonMobility(processes=processes)
            im_predictor.add_im_predictions(psm_list)

        with Encoder.from_psm_list(psm_list) as encoder:
            # Use our custom parallelized class with our spectrum reader
            custom_parallelized = PatchParallelized(
                encoder=encoder,
                model=model,
                model_dir=model_dir,
                ms2_tolerance=ms2_tolerance,
                processes=processes
            )

            logging.info("Processing spectra and peptides with custom reader...")
            results = custom_parallelized.process_spectra(
                psm_list, spectrum_file, spectrum_id_pattern
            )

            # Correlations also requested
            if compute_correlations:
                logging.info("Computing correlations")
                calculate_correlations(results)
                logging.info(f"Median correlation: {np.median(list(r.correlation for r in results))}")

            return results


def read_spectrum_file(spec_file: str) -> Generator[ObservedSpectrum, None, None]:
    """
    Read MS2 spectra from a supported file format; inferring the type from the filename extension.

    Parameters
    ----------
    spectrum_file
        Path to MGF or mzML file.

    Yields
    ------
    ObservedSpectrum

    Raises
    ------
    UnsupportedSpectrumFiletypeError
        If the file extension is not supported.

    """
    try:
        spectra = OpenMSHelper.get_mslevel_spectra(file_name=str(spec_file), ms_level=2)
    except ValueError:
        raise exceptions.UnsupportedSpectrumFiletypeError(Path(spec_file).suffixes)

    for spectrum in spectra:
        mz, intensities = spectrum.get_peaks()
        precursors = spectrum.getPrecursors()
        obs_spectrum = None
        if len(precursors) > 0:
            precursor = precursors[0]
            charge_state = precursor.getCharge()
            exp_mz = precursor.getMZ()
            rt = spectrum.getRT()
            spec_id = spectrum.getNativeID()

            obs_spectrum = ObservedSpectrum(
                mz=np.array(mz, dtype=np.float32),
                intensity=np.array(intensities, dtype=np.float32),
                identifier=str(spec_id),
                precursor_mz=float(exp_mz),
                precursor_charge=float(charge_state),
                retention_time=float(rt)
            )
        # Workaround for mobiusklein/mzdata#3
        if (
            obs_spectrum == None or
            obs_spectrum.identifier == ""
            or obs_spectrum.mz.shape[0] == 0
            or obs_spectrum.intensity.shape[0] == 0
        ):
            continue
        yield obs_spectrum


def _custom_process_spectra(
    enumerated_psm_list: List[Tuple[int, PSM]],
    spec_file: str,
    vector_file: bool,
    encoder: Encoder,
    model: str,
    ms2_tolerance: float,
    spectrum_id_pattern: str,
    annotations_only: bool = False,
) -> List[ProcessingResult]:
    """
    Perform requested tasks for each spectrum in spectrum file.

    Parameters
    ----------
    enumerated_psm_list
        List of tuples of (index, PSM) for each PSM in the input file.
    spec_file
        Filename of spectrum file
    vector_file
        If feature vectors should be extracted instead of predictions
    encoder: Encoder
        Configured encoder to use for peptide and peptidoform encoding
    model
        Name of prediction model to be used
    ms2_tolerance
        Fragmentation spectrum m/z error tolerance in Dalton
    spectrum_id_pattern
        Regular expression pattern to apply to spectrum titles before matching to
        peptide file entries
    annotations_only
        If only peak annotations should be extracted from the spectrum file

    """
    ms2pip_pyx.ms2pip_init(*encoder.encoder_files)
    results = []
    ion_types = [it.lower() for it in MODELS[model]["ion_types"]]

    try:
        spectrum_id_regex = re.compile(spectrum_id_pattern)
    except TypeError:
        spectrum_id_regex = re.compile(r"(.*)")

    # Restructure PeptideRecord entries as spec_id -> [(id, psm_1), (id, psm_2), ...]
    psms_by_specid = defaultdict(list)
    for psm_index, psm in enumerated_psm_list:
        psms_by_specid[str(psm.spectrum_id)].append((psm_index, psm))

    for spectrum in read_spectrum_file(spec_file):
        # Match spectrum ID with provided regex, use first match group as new ID
        match = spectrum_id_regex.search(spectrum.identifier)
        try:
            spectrum_id = match[1]
        except (TypeError, IndexError):
            raise exceptions.TitlePatternError(
                f"Spectrum title pattern `{spectrum_id_pattern}` could not be matched to "
                f"spectrum ID `{spectrum.identifier}`. "
                " Are you sure that the regex contains a capturing group?"
            )

        if spectrum_id not in psms_by_specid:
            continue

        # Spectrum preprocessing:
        # Remove reporter ions and precursor peak, normalize, transform
        for label_type in ["iTRAQ", "TMT"]:
            if label_type in model:
                spectrum.remove_reporter_ions(label_type)
        # spectrum.remove_precursor()  # TODO: Decide to implement this or not
        spectrum.tic_norm()
        spectrum.log2_transform()

        for psm_index, psm in psms_by_specid[spectrum_id]:
            try:
                enc_peptidoform = encoder.encode_peptidoform(psm.peptidoform)
            except exceptions.InvalidAminoAcidError:
                result = ProcessingResult(psm_index=psm_index, psm=psm)
                results.append(result)
                continue

            targets = ms2pip_pyx.get_targets(
                enc_peptidoform,
                spectrum.mz.astype(np.float32),
                spectrum.intensity.astype(np.float32),
                float(ms2_tolerance),
                MODELS[model]["peaks_version"],
            )
            targets = {i: np.array(t, dtype=np.float32) for i, t in zip(ion_types, targets)}

            if not psm.peptidoform.precursor_charge:
                psm.peptidoform.precursor_charge = spectrum.precursor_charge

            if vector_file:
                enc_peptide = encoder.encode_peptide(psm.peptidoform)
                feature_vectors = np.array(
                    ms2pip_pyx.get_vector(
                        enc_peptide, enc_peptidoform, psm.peptidoform.precursor_charge
                    ),
                    dtype=np.uint16,
                )
                result = ProcessingResult(
                    psm_index=psm_index,
                    psm=psm,
                    theoretical_mz=None,
                    predicted_intensity=None,
                    observed_intensity=targets,
                    correlation=None,
                    feature_vectors=feature_vectors,
                )

            elif annotations_only:
                # Only return mz and targets
                mz = ms2pip_pyx.get_mzs(enc_peptidoform, MODELS[model]["peaks_version"])
                mz = {i: np.array(mz, dtype=np.float32) for i, mz in zip(ion_types, mz)}

                result = ProcessingResult(
                    psm_index=psm_index,
                    psm=psm,
                    theoretical_mz=mz,
                    predicted_intensity=None,
                    observed_intensity=targets,
                    correlation=None,
                    feature_vectors=None,
                )

            else:
                # Predict with C model or get feature vectors for XGBoost
                try:
                    result = _process_peptidoform(psm_index, psm, model, encoder, ion_types)
                except (
                    exceptions.InvalidPeptidoformError,
                    exceptions.InvalidAminoAcidError,
                ):
                    result = ProcessingResult(psm_index=psm_index, psm=psm)
                else:
                    result.observed_intensity = targets

            results.append(result)

    return results



