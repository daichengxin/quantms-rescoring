import logging
from itertools import chain
from typing import Optional

from ms2pip import correlate
from ms2pip.exceptions import NoMatchingSpectraFound
from ms2rescore.feature_generators import MS2PIPFeatureGenerator
from ms2rescore.feature_generators.base import FeatureGeneratorException
from ms2rescore.utils import infer_spectrum_path
from psm_utils import PSMList

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
        lower_score_is_better: bool = True,
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
        self._calibration_set_size = calibration_set_size
        self._correlation_threshold = correlation_threshold
        self._lower_score_is_better = lower_score_is_better

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
                    lower_score_is_better=self._lower_score_is_better,
                )

                if not valid_correlation:
                    logging.error(
                        "Invalid correlation found. Please try a different model or adjust the correlation threshold."
                    )
                    raise FeatureGeneratorException(
                        f"Invalid correlation found. Please try a different model or adjust the correlation threshold."
                    )

                self._calculate_features(psm_list_run, ms2pip_results)

                current_run += 1

    def _validate_scores(
        self, ms2pip_results, calibration_set_size, correlation_threshold, lower_score_is_better
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
        lower_score_is_better : bool
            Indicates if a lower PSM score is considered better.

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
        ms2pip_results_copy.sort(key=lambda x: x.psm.score, reverse=not lower_score_is_better)

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
