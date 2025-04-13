import itertools
import multiprocessing
import re
from collections import defaultdict
from itertools import chain
from pathlib import Path
from typing import Optional, Tuple, List, Union, Callable, Generator, Dict

import ms2pip.exceptions as exceptions
import numpy as np
from math import ceil
from ms2pip._cython_modules import ms2pip_pyx
from ms2pip._utils.encoder import Encoder
from ms2pip._utils.ion_mobility import IonMobility
from ms2pip._utils.psm_input import read_psms
from ms2pip._utils.retention_time import RetentionTime
from ms2pip.constants import MODELS
from ms2pip.core import _Parallelized, _process_peptidoform
from ms2pip.exceptions import NoMatchingSpectraFound
from ms2pip.result import ProcessingResult
from ms2pip.spectrum import ObservedSpectrum
from ms2rescore.feature_generators import MS2PIPFeatureGenerator
from ms2rescore.feature_generators.base import FeatureGeneratorException
from ms2rescore.utils import infer_spectrum_path
from psm_utils import PSMList, PSM

from quantmsrescore.constants import SUPPORTED_MODELS_MS2PIP
from quantmsrescore.exceptions import Ms2pipIncorrectModelException
from quantmsrescore.logging_config import get_logger
from quantmsrescore.openms import OpenMSHelper

# Get logger for this module
logger = get_logger(__name__)


class PatchParallelized(_Parallelized):
    """
    Extended version of _Parallelized that supports custom spectrum file reading using pyopenms instead of
    ms2rescore_rs. We need to read the files using OpenMS mzML handler rather than mzdata tools in rust because if
    the mzML doesn't work on this component, it will not work on quantms but if it not works on mzdata is not in our control
    to changed and those fails may be supported and OK in quantms.
    See example in issue:
    """

    def __init__(self, encoder, model=None, model_dir=None, ms2_tolerance=0.02, processes=None):
        """
        Initialize with all original parameters plus a custom spectrum reader.

        Parameters
        ----------
        encoder:
            Encoder object to use for encoding peptides.
        model:
            MS²PIP model to use for scoring.
        model_dir:
            Directory containing model files.
        ms2_tolerance:
            MS² tolerance in Da.
        processes:
            Number of processes to use for parallelization.
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
        logger.debug(f"Gathered data for {len(results)} PSMs.")

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
                logger.warning("No PSMs to process.")
                return []

            logger.info(
                "The Pool number of process {} and CPUs {}".format(
                    pool._processes, multiprocessing.cpu_count()
                )
            )

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

            logger.info(f"Processing {len(chunks)} chunk(s) of ~{chunk_size} entries each.")

            # Add jobs to pool
            mp_results = [
                pool.apply_async(func, args=(psm_list_chunk, *args)) for psm_list_chunk in chunks
            ]
            results = [r.get() for r in mp_results]

            pool.close()
            pool.join()

        # Sort results by input order
        results = sorted(
            itertools.chain.from_iterable(results), key=lambda result: result.psm_index
        )

        return results

    def _get_pool(self):
        """Get multiprocessing pool with recursion protection."""
        logger.debug(f"Starting workers (processes={self.processes})...")

        if multiprocessing.current_process().daemon:
            logger.warning(
                "Running in a daemon process. Disabling multiprocessing as daemonic processes cannot have children."
            )
            return multiprocessing.dummy.Pool(1)

        if self.processes == 1:
            logger.debug("Using dummy multiprocessing pool.")
            return multiprocessing.dummy.Pool(1)

        # Check if already inside a worker process
        if multiprocessing.parent_process() is not None:
            logger.warning(
                "Attempting to create a pool inside a worker process! Returning a dummy pool instead."
            )
            return multiprocessing.dummy.Pool(1)

        return multiprocessing.get_context("spawn").Pool(self.processes)


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
        valid_correlations_size: Optional[float] = 0.70,
        correlation_threshold: Optional[float] = 0.6,
        higher_score_better: bool = True,
        force_model: bool = False,
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
        self._calibration_set_size: float = calibration_set_size
        self._valid_correlations_size: float = valid_correlations_size
        self._correlation_threshold: float = correlation_threshold
        self._higher_score_better: bool = higher_score_better
        self._force_model: bool = force_model

    def add_features(self, psm_list: PSMList) -> None:
        """
        Add MS²PIP-derived features to PSMs.

        Parameters
        ----------
        psm_list
        PSMs to add features to.
        """
        logger.info("Adding MS²PIP-derived features to PSMs.")
        psm_dict = psm_list.get_psm_dict()
        current_run = 1
        total_runs = sum(len(runs) for runs in psm_dict.values())

        for runs in psm_dict.values():
            for run, psms in runs.items():
                logger.info(
                    f"Running MS²PIP {self.model} for PSMs from run ({current_run}/{total_runs}) `{run}`..."
                )
                psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
                spectrum_filename = infer_spectrum_path(self.spectrum_path, run)
                logger.debug(f"Using spectrum file `{spectrum_filename}`")
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
                    logger.error(
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

    def _validate_scores(
        self,
        ms2pip_results,
        calibration_set_size,
        valid_correlations_size,
        correlation_threshold,
        higher_score_better,
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
        ms2pip_results_copy = [
            result
            for result in ms2pip_results_copy
            if not result.psm.is_decoy and result.psm.rank == 1
        ]
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

        logger.info(
            f"The percentage of PSMs in the top {calibration_set_size * 100}% with a correlation greater than {correlation_threshold} is: "
            f"{(len(valid_correlation) / len(calibration_set)) * 100:.2f}%"
        )

        if not self._force_model:
            # If the number of valid correlations is less than 80% of the calibration set, return False
            if len(valid_correlation) < len(calibration_set) * valid_correlations_size:
                return False

        return True

    def _find_best_ms2pip_model(
        self, batch_psms: PSMList, knwon_fragmentation: Optional[str] = None
    ) -> Tuple[str, float]:
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

        for fragment_types in filtered_models:
            for model in filtered_models[fragment_types]:
                logger.info(f"Running MS²PIP for model `{model}`...")
                ms2pip_results = self.custom_correlate(
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

        return best_model, best_correlation

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
        total_correlation = sum(
            [
                psm.correlation
                for psm in ms2pip_results
                if psm.correlation is not None and not np.isnan(psm.correlation)
            ]
        )
        return total_correlation / len(ms2pip_results)

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
            logger.info("Adding retention time predictions")
            rt_predictor = RetentionTime(processes=processes)
            rt_predictor.add_rt_predictions(psm_list)

        if add_ion_mobility:
            logger.info("Adding ion mobility predictions")
            im_predictor = IonMobility(processes=processes)
            im_predictor.add_im_predictions(psm_list)

        with Encoder.from_psm_list(psm_list) as encoder:
            # Use our custom parallelized class with our spectrum reader
            custom_parallelized = PatchParallelized(
                encoder=encoder,
                model=model,
                model_dir=model_dir,
                ms2_tolerance=ms2_tolerance,
                processes=processes,
            )

            logger.info("Processing spectra and peptides with custom reader...")
            results = custom_parallelized.process_spectra(
                psm_list, spectrum_file, spectrum_id_pattern
            )

            # Correlations also requested
            if compute_correlations:
                calculate_correlations(results)
                logger.info(
                    f"Median correlation: {np.median([r.correlation for r in results if r.correlation is not None and not np.isnan(r.correlation)])}, model {model}"
                )

            return results


def calculate_correlations(results: List[ProcessingResult]) -> None:
    """Calculate and add Pearson correlations to list of results."""
    for result in results:
        if result.predicted_intensity and result.observed_intensity:
            pred_int = np.concatenate([i for i in result.predicted_intensity.values()])
            obs_int = np.concatenate([i for i in result.observed_intensity.values()])
            result.correlation = np.corrcoef(pred_int, obs_int)[0][1]
        else:
            result.correlation = None
            logger.info("Results {} is empty".format(result))


def read_spectrum_file(spec_file: str) -> Generator[ObservedSpectrum, None, None]:
    """
    Read MS2 spectra from a supported file format; inferring the type from the filename extension.

    Parameters
    ----------
    spec_file:
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
                retention_time=float(rt),
            )
        if (
            obs_spectrum is None
            or obs_spectrum.identifier == ""
            or obs_spectrum.mz.shape[0] == 0
            or obs_spectrum.intensity.shape[0] == 0
        ):
            continue
        yield obs_spectrum


def _organize_psms_by_spectrum_id(
    enumerated_psm_list: List[Tuple[int, PSM]]
) -> Dict[str, List[Tuple[int, PSM]]]:
    """
    Organize PSMs by spectrum ID for efficient lookup.

    Parameters
    ----------
    enumerated_psm_list
        List of tuples of (index, PSM) for each PSM in the input file.

    Returns
    -------
    Dict[str, List[Tuple[int, PSM]]]
        Dictionary mapping spectrum IDs to lists of (index, PSM) tuples.
    """
    psms_by_specid = defaultdict(list)
    for psm_index, psm in enumerated_psm_list:
        psms_by_specid[str(psm.spectrum_id)].append((psm_index, psm))
    return psms_by_specid


def _preprocess_spectrum(spectrum: ObservedSpectrum, model: str) -> None:
    """
    Preprocess a spectrum by removing reporter ions, normalizing, and transforming.

    Parameters
    ----------
    spectrum
        The spectrum to preprocess.
    model
        The model name, used to determine if reporter ions should be removed.
    """
    # Remove reporter ions if needed
    for label_type in ["iTRAQ", "TMT"]:
        if label_type in model:
            spectrum.remove_reporter_ions(label_type)

    spectrum.tic_norm()
    spectrum.log2_transform()


def _get_targets_for_psm(
    psm: PSM,
    spectrum: ObservedSpectrum,
    encoder: Encoder,
    ms2_tolerance: float,
    model: str,
    ion_types: List[str],
) -> Tuple[Optional[np.ndarray], Dict[str, np.ndarray]]:
    """
    Get targets for a PSM from a spectrum.

    Parameters
    ----------
    psm
        The PSM to get targets for.
    spectrum
        The spectrum to get targets from.
    encoder
        The encoder to use for peptide and peptidoform encoding.
    ms2_tolerance
        The MS2 tolerance to use.
    model
        The model name.
    ion_types
        The ion types to use.

    Returns
    -------
    Tuple[Optional[np.ndarray], Dict[str, np.ndarray]]
        A tuple containing the encoded peptidoform and the targets.
    """
    try:
        enc_peptidoform = encoder.encode_peptidoform(psm.peptidoform)
    except exceptions.InvalidAminoAcidError:
        return None, {}

    # Get targets
    targets = ms2pip_pyx.get_targets(
        enc_peptidoform,
        spectrum.mz.astype(np.float32),
        spectrum.intensity.astype(np.float32),
        float(ms2_tolerance),
        MODELS[model]["peaks_version"],
    )
    targets = {i: np.array(t, dtype=np.float32) for i, t in zip(ion_types, targets)}

    # Set precursor charge if not set
    if not psm.peptidoform.precursor_charge:
        psm.peptidoform.precursor_charge = spectrum.precursor_charge

    return enc_peptidoform, targets


def _create_result_for_mode(
    psm_index: int,
    psm: PSM,
    enc_peptidoform: np.ndarray,
    targets: Dict[str, np.ndarray],
    vector_file: bool,
    annotations_only: bool,
    model: str,
    encoder: Encoder,
    ion_types: List[str],
) -> ProcessingResult:
    """
    Create a ProcessingResult based on the processing mode.

    Parameters
    ----------
    psm_index
        The index of the PSM.
    psm
        The PSM.
    enc_peptidoform
        The encoded peptidoform.
    targets
        The targets.
    vector_file
        Whether to extract feature vectors.
    annotations_only
        Whether to extract only annotations.
    model
        The model name.
    encoder
        The encoder.
    ion_types
        The ion types.

    Returns
    -------
    ProcessingResult
        The processing result.
    """
    if vector_file:
        # Extract feature vectors
        enc_peptide = encoder.encode_peptide(psm.peptidoform)
        feature_vectors = np.array(
            ms2pip_pyx.get_vector(enc_peptide, enc_peptidoform, psm.peptidoform.precursor_charge),
            dtype=np.uint16,
        )
        return ProcessingResult(
            psm_index=psm_index,
            psm=psm,
            theoretical_mz=None,
            predicted_intensity=None,
            observed_intensity=targets,
            correlation=None,
            feature_vectors=feature_vectors,
        )
    elif annotations_only:
        # Extract only annotations
        mz = ms2pip_pyx.get_mzs(enc_peptidoform, MODELS[model]["peaks_version"])
        mz = {i: np.array(mz, dtype=np.float32) for i, mz in zip(ion_types, mz)}

        return ProcessingResult(
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

        return result


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

    Returns
    -------
    List[ProcessingResult]
        List of processing results.
    """
    # Initialize MS2PIP
    ms2pip_pyx.ms2pip_init(*encoder.encoder_files)
    results = []
    ion_types = [it.lower() for it in MODELS[model]["ion_types"]]

    # Compile regex for spectrum ID matching
    try:
        spectrum_id_regex = re.compile(spectrum_id_pattern)
    except TypeError:
        spectrum_id_regex = re.compile(r"(.*)")

    # Organize PSMs by spectrum ID
    psms_by_specid = _organize_psms_by_spectrum_id(enumerated_psm_list)

    # Process each spectrum
    for spectrum in read_spectrum_file(spec_file):
        # Match spectrum ID with provided regex
        match = spectrum_id_regex.search(spectrum.identifier)
        try:
            spectrum_id = match[1]
        except (TypeError, IndexError):
            raise exceptions.TitlePatternError(
                f"Spectrum title pattern `{spectrum_id_pattern}` could not be matched to "
                f"spectrum ID `{spectrum.identifier}`. "
                " Are you sure that the regex contains a capturing group?"
            )

        # Skip if no matching PSMs
        if spectrum_id not in psms_by_specid:
            continue

        # Preprocess spectrum
        _preprocess_spectrum(spectrum, model)

        # Process each PSM for this spectrum
        for psm_index, psm in psms_by_specid[spectrum_id]:
            # Get targets for PSM
            enc_peptidoform, targets = _get_targets_for_psm(
                psm, spectrum, encoder, ms2_tolerance, model, ion_types
            )

            # Skip if encoding failed
            if enc_peptidoform is None:
                result = ProcessingResult(psm_index=psm_index, psm=psm)
                results.append(result)
                continue

            # Create result based on processing mode
            result = _create_result_for_mode(
                psm_index,
                psm,
                enc_peptidoform,
                targets,
                vector_file,
                annotations_only,
                model,
                encoder,
                ion_types,
            )

            results.append(result)

    return results
