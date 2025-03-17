import logging
from itertools import chain
from typing import Optional, Tuple, List

from ms2pip import correlate
from ms2pip.exceptions import NoMatchingSpectraFound
from ms2pip.result import ProcessingResult
from ms2rescore.feature_generators import MS2PIPFeatureGenerator
from ms2rescore.feature_generators.base import FeatureGeneratorException
from ms2rescore.utils import infer_spectrum_path
from psm_utils import PSMList

from quantmsrescore.constants import SUPPORTED_MODELS_MS2PIP
from quantmsrescore.exceptions import Ms2pipIncorrectModelException

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


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
                    f"Running MS²PIP for PSMs from run ({current_run}/{total_runs}) `{run}`..."
                )
                psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))
                spectrum_filename = infer_spectrum_path(self.spectrum_path, run)
                logging.debug(f"Using spectrum file `{spectrum_filename}`")
                try:
                    ms2pip_results = correlate(
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
                    correlation_threshold=self._correlation_threshold,
                    higher_score_better=self._higher_score_better,
                )

                if not valid_correlation:
                    logging.error(
                        "Invalid correlation found. Please try a different model or adjust the correlation threshold."
                    )
                    raise Ms2pipIncorrectModelException(
                        message="Invalid correlation found. Please try a different model or adjust the correlation threshold.",
                        model=self.model,
                    )
                self._calculate_features(psm_list_run, ms2pip_results)
                current_run += 1

    @staticmethod
    def _validate_scores(
        ms2pip_results, calibration_set_size, correlation_threshold, higher_score_better
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

        # If the number of valid correlations is less than 80% of the calibration set, return False
        if len(valid_correlation) < len(calibration_set) * 0.8:
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
