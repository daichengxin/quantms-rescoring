import logging
from itertools import chain

from ms2pip import correlate
from ms2pip.exceptions import NoMatchingSpectraFound
from ms2rescore.feature_generators import MS2PIPFeatureGenerator
from ms2rescore.feature_generators.base import FeatureGeneratorException
from ms2rescore.utils import infer_spectrum_path
from psm_utils import PSMList

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


class MS2PIPAnnotator(MS2PIPFeatureGenerator):

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
                self._calculate_features(psm_list_run, ms2pip_results)
                current_run += 1
