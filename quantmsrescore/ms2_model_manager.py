import pandas as pd
from peptdeep.pretrained_models import ModelManager, _download_models, model_mgr_settings, MODEL_ZIP_FILE_PATH, \
    psm_sampling_with_important_mods
from peptdeep.model.ms2 import pDeepModel
from peptdeep.model.rt import AlphaRTModel
from peptdeep.model.ccs import AlphaCCSModel
from peptdeep.model.charge import ChargeModelForModAASeq
import os
from peptdeep.utils import logging


class MS2ModelManager(ModelManager):
    def __init__(self,
                 mask_modloss: bool = False,
                 device: str = "gpu",
                 model_dir: str = None,
                 ):
        self._train_psm_logging = True

        self.ms2_model: pDeepModel = pDeepModel(
            mask_modloss=mask_modloss, device=device
        )
        self.rt_model: AlphaRTModel = AlphaRTModel(device=device)
        self.ccs_model: AlphaCCSModel = AlphaCCSModel(device=device)

        self.charge_model: ChargeModelForModAASeq = ChargeModelForModAASeq(
            device=device
        )

        if model_dir is not None and os.path.exists(os.path.join(model_dir, "ms2.pth")):
            self.load_external_models(ms2_model_file=os.path.join(model_dir, "ms2.pth"))
            self.model_str = model_dir
        else:
            _download_models(MODEL_ZIP_FILE_PATH)
            self.load_installed_models()
            self.model_str = "generic"
        self.reset_by_global_settings(reload_models=False)

    def __str__(self):
        return self.model_str

    def ms2_fine_tuning(self, psms_df: pd.DataFrame,
                        match_intensity_df: pd.DataFrame,
                        psm_num_to_train_ms2: int = 100000000,
                        use_grid_nce_search: bool = False,
                        top_n_mods_to_train: int = 10,
                        psm_num_per_mod_to_train_ms2: int = 50,
                        psm_num_to_test_ms2: int = 0,
                        epoch_to_train_ms2: int = 20,
                        train_verbose: bool = False):

        self.psm_num_to_train_ms2 = psm_num_to_train_ms2
        self.use_grid_nce_search = use_grid_nce_search
        self.top_n_mods_to_train = top_n_mods_to_train
        self.psm_num_per_mod_to_train_ms2 = psm_num_per_mod_to_train_ms2
        self.psm_num_to_test_ms2 = psm_num_to_test_ms2
        self.train_verbose = train_verbose
        self.epoch_to_train_ms2 = epoch_to_train_ms2
        self.train_ms2_model(psms_df, match_intensity_df)

    def train_ms2_model(
            self,
            psm_df: pd.DataFrame,
            matched_intensity_df: pd.DataFrame,
    ):
        """
        Using matched_intensity_df to train/fine-tune the ms2 model.

        1. It will sample `n=self.psm_num_to_train_ms2` PSMs into training dataframe (`tr_df`) to for fine-tuning.
        2. This method will also consider some important PTMs (`n=self.top_n_mods_to_train`) into `tr_df` for fine-tuning.
        3. If `self.use_grid_nce_search==True`, this method will call `self.ms2_model.grid_nce_search` to find the best NCE and instrument.

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM dataframe for fine-tuning

        matched_intensity_df : pd.DataFrame
            The matched fragment intensities for `psm_df`.
        """
        if self.psm_num_to_train_ms2 > 0:
            if self.psm_num_to_train_ms2 < len(psm_df):
                tr_df = psm_sampling_with_important_mods(
                    psm_df,
                    self.psm_num_to_train_ms2,
                    self.top_n_mods_to_train,
                    self.psm_num_per_mod_to_train_ms2,
                ).copy()
            else:
                tr_df = psm_df
            if len(tr_df) > 0:
                tr_inten_df = pd.DataFrame()
                for frag_type in self.ms2_model.charged_frag_types:
                    if frag_type in matched_intensity_df.columns:
                        tr_inten_df[frag_type] = matched_intensity_df[frag_type]
                    else:
                        tr_inten_df[frag_type] = 0.0
                # normalize_fragment_intensities(tr_df, tr_inten_df)

                if self.use_grid_nce_search:
                    self.nce, self.instrument = self.ms2_model.grid_nce_search(
                        tr_df,
                        tr_inten_df,
                        nce_first=model_mgr_settings["transfer"]["grid_nce_first"],
                        nce_last=model_mgr_settings["transfer"]["grid_nce_last"],
                        nce_step=model_mgr_settings["transfer"]["grid_nce_step"],
                        search_instruments=model_mgr_settings["transfer"][
                            "grid_instrument"
                        ],
                    )
                    tr_df["nce"] = self.nce
                    tr_df["instrument"] = self.instrument
                else:
                    self.set_default_nce_instrument(tr_df)
        else:
            tr_df = []

        if self.psm_num_to_test_ms2 > 0:
            if len(tr_df) > 0:
                test_psm_df = psm_df[~psm_df.sequence.isin(set(tr_df.sequence))].copy()
                if len(test_psm_df) > self.psm_num_to_test_ms2:
                    test_psm_df = test_psm_df.sample(n=self.psm_num_to_test_ms2)
                elif len(test_psm_df) == 0:
                    logging.info(
                        "No enough PSMs for testing MS2 models, "
                        "please reduce the `psm_num_to_train_ms2` "
                        "value according to overall PSM numbers. "
                    )
                    test_psm_df = []
            else:
                test_psm_df = psm_df.copy()
                tr_inten_df = pd.DataFrame()
                for frag_type in self.ms2_model.charged_frag_types:
                    if frag_type in matched_intensity_df.columns:
                        tr_inten_df[frag_type] = matched_intensity_df[frag_type]
                    else:
                        tr_inten_df[frag_type] = 0.0
            self.set_default_nce_instrument(test_psm_df)
        else:
            test_psm_df = []

        if len(test_psm_df) > 0:
            logging.info(
                "Testing pretrained MS2 model on testing df:\n"
                + str(self.ms2_model.test(test_psm_df, tr_inten_df))
            )
        if len(tr_df) > 0:
            if self._train_psm_logging:
                logging.info(
                    f"{len(tr_df)} PSMs for MS2 model training/transfer learning"
                )
            self.ms2_model.train(
                tr_df,
                fragment_intensity_df=tr_inten_df,
                batch_size=self.batch_size_to_train_ms2,
                epoch=self.epoch_to_train_ms2,
                warmup_epoch=self.warmup_epoch_to_train_ms2,
                lr=self.lr_to_train_ms2,
                verbose=self.train_verbose,
            )
            logging.info(
                "Testing refined MS2 model on training df:\n"
                + str(self.ms2_model.test(tr_df, tr_inten_df))
            )
        if len(test_psm_df) > 0:
            logging.info(
                "Testing refined MS2 model on testing df:\n"
                + str(self.ms2_model.test(test_psm_df, tr_inten_df))
            )

        self.model_str = "retrained_model"

    def save_ms2_model(self):
        self.ms2_model.save("retrained_ms2.pth")
