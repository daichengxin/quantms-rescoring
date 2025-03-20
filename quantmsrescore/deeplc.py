import contextlib
import os
from collections import defaultdict
from itertools import chain

import numpy as np
from ms2rescore.feature_generators import DeepLCFeatureGenerator
from psm_utils import PSMList

from quantmsrescore.logging_config import get_logger

# Get logger for this module
logger = get_logger(__name__)


class DeepLCAnnotator(DeepLCFeatureGenerator):

    def add_features(self, psm_list: PSMList) -> None:
        """Add DeepLC-derived features to PSMs."""

        logger.info("Adding DeepLC-derived features to PSMs.")

        # Get easy-access nested version of PSMList
        psm_dict = psm_list.get_psm_dict()

        # Run DeepLC for each spectrum file
        current_run = 1
        total_runs = sum(len(runs) for runs in psm_dict.values())

        for runs in psm_dict.values():
            # Reset DeepLC predictor for each collection of runs
            self.deeplc_predictor = None
            self.selected_model = None
            for run, psms in runs.items():
                peptide_rt_diff_dict = defaultdict(
                    lambda: {
                        "observed_retention_time_best": np.inf,
                        "predicted_retention_time_best": np.inf,
                        "rt_diff_best": np.inf,
                    }
                )
                logger.info(
                    f"Running DeepLC for PSMs from run ({current_run}/{total_runs}): `{run}`..."
                )

                # Disable wild logging to stdout by Tensorflow, unless in debug mode
                with (
                    contextlib.redirect_stdout(open(os.devnull, "w", encoding="utf-8"))
                    if not self._verbose and os.devnull is not None
                    else contextlib.nullcontext()
                ):
                    # Make new PSM list for this run (chain PSMs per spectrum to flat list)
                    psm_list_run = PSMList(psm_list=list(chain.from_iterable(psms.values())))

                    psm_list_calibration = self._get_calibration_psms(psm_list_run)
                    logger.debug(f"Calibrating DeepLC with {len(psm_list_calibration)} PSMs...")
                    self.deeplc_predictor = self.DeepLC(
                        n_jobs=self.processes,
                        verbose=self._verbose,
                        path_model=self.selected_model or self.user_model,
                        **self.deeplc_kwargs,
                    )
                    self.deeplc_predictor.calibrate_preds(
                        psm_list=psm_list_calibration, return_plotly_report=True
                    )
                    # Still calibrate for each run, but do not try out all model options.
                    # Just use the model that was selected based on the first run
                    self.selected_model = list(self.deeplc_predictor.model.keys())
                    logger.debug(
                        f"Selected DeepLC model {self.selected_model} based on "
                        "calibration of first run. Using this model (after new "
                        "calibrations) for the remaining runs."
                    )
                    logger.debug("Predicting retention times...")
                    predictions = np.array(self.deeplc_predictor.make_preds(psm_list_run))
                    observations = psm_list_run["retention_time"]
                    rt_diffs_run = np.abs(predictions - observations)

                    logger.debug("Adding features to PSMs...")
                    for i, psm in enumerate(psm_list_run):
                        psm["rescoring_features"].update(
                            {
                                "observed_retention_time": observations[i],
                                "predicted_retention_time": predictions[i],
                                "rt_diff": rt_diffs_run[i],
                            }
                        )
                        peptide = psm.peptidoform.proforma.split("\\")[0]  # remove charge
                        if peptide_rt_diff_dict[peptide]["rt_diff_best"] > rt_diffs_run[i]:
                            peptide_rt_diff_dict[peptide] = {
                                "observed_retention_time_best": observations[i],
                                "predicted_retention_time_best": predictions[i],
                                "rt_diff_best": rt_diffs_run[i],
                            }
                    for psm in psm_list_run:
                        psm["rescoring_features"].update(
                            peptide_rt_diff_dict[psm.peptidoform.proforma.split("\\")[0]]
                        )
                current_run += 1
