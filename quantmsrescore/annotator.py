import copy
import logging
from pathlib import Path
from typing import Union

from psm_utils import PSMList

from quantmsrescore.deeplc import DeepLCAnnotator
from quantmsrescore.exceptions import Ms2pipIncorrectModelException
from quantmsrescore.idxmlreader import IdXMLRescoringReader
from quantmsrescore.ms2pip import MS2PIPAnnotator
from quantmsrescore.openms import OpenMSHelper

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


class Annotator:
    def __init__(
        self,
        feature_generators: str,
        only_features: str = None,
        ms2pip_model: str = "HCD2021",
        ms2pip_model_path: str = "models",
        ms2_tolerance: float = 0.05,
        calibration_set_size: float = 0.2,
        skip_deeplc_retrain: bool = False,
        processes: int = 2,
        id_decoy_pattern: str = "^DECOY_",
        lower_score_is_better: bool = True,
        log_level: str = "INFO",
        spectrum_id_pattern: str = "(.*)",  # default for openms idXML
        psm_id_pattern: str = "(.*)",  # default for openms idXML
        remove_missing_spectra: bool = True,
        ms2_only: bool = True,
        find_best_ms2pip_model: bool = False,
    ):
        """
        Initializes the Annotator class with configuration parameters for feature generation
        and rescoring of peptide-spectrum matches (PSMs).

        Parameters
        ----------
        feature_generators : str
            Comma-separated list of feature generators to use for annotation (e.g., "ms2pip,deeplc").
        only_features : str, optional
            Comma-separated list of features to include in the annotation (default: None).
        ms2pip_model : str, optional
            The MS2PIP model to use for feature generation (default: "HCD2021").
        ms2pip_model_path : str, optional
            The path to the MS2PIP model directory (default: "models").
        ms2_tolerance : float, optional
            The MS2 tolerance to use for feature generation (default: 0.05).
        calibration_set_size : float, optional
            The percentage of PSMs to use for calibration (default: 0.2).
        skip_deeplc_retrain : bool, optional
            Whether to skip retraining the deepLC model (default: False).
        processes : int, optional
            The number of processes to use for feature generation (default: 2).
        id_decoy_pattern : str, optional
            The pattern to use for identifying decoy PSMs (default: "^DECOY_").
        lower_score_is_better : bool, optional
            Whether a lower score is better for feature generation (default: True).
        log_level : str, optional
            The logging level to use for the Annotator (default: "INFO").
        spectrum_id_pattern : str, optional
            The pattern to use for identifying spectrum IDs (default: "(.*)").
        psm_id_pattern : str, optional
            The pattern to use for identifying PSM IDs (default: "(.*)").
        remove_missing_spectra : bool, optional
            Whether to remove PSMs with missing spectra (default: True).
        ms2_only : bool, optional
            Whether to only process MS2-level PSMs (default: True).
        find_best_ms2pip_model : bool, optional
            Whether to find the best MS2PIP model based on the top calibration set (default: False).

        Raises
        -------
        ValueError: If no feature generators are provided.
        """
        self._idxml_reader = None
        if not feature_generators:
            raise ValueError("feature_annotators must be provided.")
        if "deeplc" not in feature_generators and "ms2pip" not in feature_generators:
            raise ValueError("At least one of deeplc or ms2pip must be provided.")

        feature_annotators = feature_generators.split(",")
        if "deeplc" in feature_annotators:
            self._deepLC = True
        else:
            self._deepLC = False
        if "ms2pip" in feature_annotators:
            self._ms2pip = True
        else:
            self._ms2pip = False

        self._only_features = []
        if only_features is not None:
            self._only_features = OpenMSHelper.validate_features(only_features.split(","))

        self._ms2pip_model = ms2pip_model
        self._ms2pip_model_path = ms2pip_model_path
        self._ms2_tolerance = ms2_tolerance
        self._calibration_set_size = calibration_set_size
        self._processes = processes
        self._id_decoy_pattern = id_decoy_pattern
        self._lower_score_is_better = lower_score_is_better
        self._log_level = log_level
        self._spectrum_id_pattern = spectrum_id_pattern
        self._psm_id_pattern = psm_id_pattern
        self._skip_deeplc_retrain = skip_deeplc_retrain
        self._remove_missing_spectra = remove_missing_spectra
        self._ms2_only = ms2_only
        self._find_best_ms2pip_model = find_best_ms2pip_model

    def build_idxml_data(self, idxml_file: Union[str, Path], spectrum_path: Union[str, Path]):
        """
        Build and load deeplc_models from an idXML file for annotation.

        This method initializes the IdXMLRescoringReader with the specified
        idXML and mzML files, and counts the number of decoy and target PSMs.

        Parameters
        ----------
        idxml_file : Union[str, Path]
            The path to the idXML file to be processed.
        spectrum_path : Union[str, Path]
            The path to the corresponding mzML file.
        """

        logging.info("Running the Annotator on file: %s", idxml_file)

        openms_helper = OpenMSHelper()

        try:
            # Load the idXML file and the corresponding mzML file
            self._idxml_reader = IdXMLRescoringReader(
                idexml_filename=idxml_file,
                mzml_file=spectrum_path,
                only_ms2=self._ms2_only,
                remove_missing_spectrum=self._remove_missing_spectra,
            )
        except Exception as e:
            logging.error(f"Failed to load input files: {str(e)}")
            raise

        psm_list = self._idxml_reader.psms
        decoys, targets = openms_helper.count_decoys_targets(self._idxml_reader.oms_peptides)
        logging.info(
            "Loaded %s PSMs from %s, %s decoys and %s targets",
            len(psm_list),
            idxml_file,
            decoys,
            targets,
        )

    def annotate(self):
        if self._idxml_reader is None:
            logging.error("No idXML deeplc_models loaded. Call build_idxml_data() first.")
            raise ValueError("No idXML deeplc_models loaded. Call build_idxml_data() first.")

        logging.debug(f"Running Annotations with following configurations: {self.__dict__}")

        if self._ms2pip:
            logging.info("Running MS2PIP on the PSMs")

            try:
                ms2pip_generator = MS2PIPAnnotator(
                    ms2_tolerance=self._ms2_tolerance,
                    model=self._ms2pip_model,
                    spectrum_path=self._idxml_reader.spectrum_path,
                    spectrum_id_pattern=self._spectrum_id_pattern,
                    model_dir=self._ms2pip_model_path,
                    calibration_set_size=self._calibration_set_size,
                    correlation_threshold=0.7,
                    lower_score_is_better=self._lower_score_is_better,
                    processes=self._processes,
                )
            except Exception as e:
                logging.error(f"Failed to initialize MS2PIP: {str(e)}")
                raise

            psm_list = self._idxml_reader.psms
            try:
                ms2pip_generator.add_features(psm_list)
                self._idxml_reader.psms = psm_list
                logging.info("MS2PIP Annotations added to the PSMs")
            except Ms2pipIncorrectModelException:
                if self._find_best_ms2pip_model:
                    logging.info(
                        "Finding best MS2PIP model - for now is a brute force search on the top calibrarion set"
                    )
                    batch_psms = self._get_top_batch_psms(psm_list)
                    model, corr = ms2pip_generator._find_best_ms2pip_model(
                        batch_psms=batch_psms,
                        knwon_fragmentation=self._get_highest_fragmentation(),
                    )
                    if model is not None:
                        logging.info(f"Best model found: {model} with average correlation {corr}")
                        ms2pip_generator = MS2PIPAnnotator(
                            ms2_tolerance=self._ms2_tolerance,
                            model=model,
                            spectrum_path=self._idxml_reader.spectrum_path,
                            spectrum_id_pattern=self._spectrum_id_pattern,
                            model_dir=self._ms2pip_model_path,
                            calibration_set_size=self._calibration_set_size,
                            correlation_threshold=0.7,
                            lower_score_is_better=self._lower_score_is_better,
                            processes=self._processes,
                        )
                        ms2pip_generator.add_features(psm_list)
                        self._idxml_reader.psms = psm_list
                        logging.info("MS2PIP Annotations added to the PSMs")
                    else:
                        logging.error(
                            "Not good model found for this deeplc_models please review parameters"
                        )
            except Exception as e:
                logging.error(f"Failed to add MS2PIP features: {str(e)}")

        if self._deepLC:
            logging.info("Running deepLC on the PSMs")

            try:
                kwargs = {}
                if self._skip_deeplc_retrain:
                    deeplc_annotator = DeepLCAnnotator(
                        self._lower_score_is_better,
                        calibration_set_size=self._calibration_set_size,
                        processes=self._processes,
                        **kwargs,
                    )
                else:
                    kwargs = {
                        "deeplc_retrain": True,
                    }
                    retrain_deeplc_annotator = DeepLCAnnotator(
                        self._lower_score_is_better,
                        calibration_set_size=self._calibration_set_size,
                        processes=self._processes,
                        **kwargs,
                    )
                    temp_psms = self._idxml_reader.psms.psm_list
                    retrained_psms = PSMList(psm_list=copy.deepcopy(temp_psms))
                    retrain_deeplc_annotator.add_features(retrained_psms)

                    kwargs = {
                        "deeplc_retrain": False,
                    }
                    deeplc_annotator_model = DeepLCAnnotator(
                        self._lower_score_is_better,
                        calibration_set_size=self._calibration_set_size,
                        processes=self._processes,
                        **kwargs,
                    )
                    psms = PSMList(psm_list=copy.deepcopy(temp_psms))
                    deeplc_annotator_model.add_features(psms)

                    # Compare the results of the retrained and non-retrained models
                    mae_trained = self._get_mae_from_psm_list(retrained_psms)
                    mae_model = self._get_mae_from_psm_list(psms)

                    if mae_trained < mae_model:
                        logging.info(
                            f"Retrained DeepLC model has lower MAE, using it. {retrain_deeplc_annotator.selected_model}"
                        )
                        deeplc_annotator = retrain_deeplc_annotator
                    else:
                        logging.info(
                            f"Retrained DeepLC model has higher MAE, using the original model. {deeplc_annotator_model.selected_model}"
                        )
                        deeplc_annotator = deeplc_annotator_model

            except Exception as e:
                logging.error(f"Failed to initialize DeepLC: {str(e)}")
                raise

            psm_list = self._idxml_reader.psms
            deeplc_annotator.add_features(psm_list)
            self._idxml_reader.psms = psm_list

        if self._ms2pip or self._deepLC:
            self._convert_features_psms_to_oms_peptides()

        logging.info("Annotations added to the PSMs, starting to modified OMS peptides")

    def write_idxml_file(self, filename: Union[str, Path]):
        try:
            OpenMSHelper.write_idxml_file(
                filename=filename,
                protein_ids=self._idxml_reader.openms_proteins,
                peptide_ids=self._idxml_reader.openms_peptides,
            )
            logging.info("Annotated idXML file written to %s", filename)
        except Exception as e:
            logging.error(f"Failed to write annotated idXML file: {str(e)}")
            raise

    def _convert_features_psms_to_oms_peptides(self):

        psm_dict = {next(iter(psm.provenance_data)): psm for psm in self._idxml_reader.psms}

        oms_peptides = []

        features = []
        for oms_peptide in self._idxml_reader.oms_peptides:
            hits = []
            for oms_psm in oms_peptide.getHits():
                psm_hash = OpenMSHelper.get_psm_hash_unique_id(
                    peptide_hit=oms_peptide, psm_hit=oms_psm
                )
                psm = psm_dict.get(psm_hash)
                if psm is None:
                    logging.warning(f"PSM not found for peptide {oms_peptide.getMetaValue('id')}")
                else:
                    for feature, value in psm.rescoring_features.items():
                        canonical_feature = OpenMSHelper.get_canonical_feature(feature)
                        if canonical_feature is not None:
                            if (
                                self._only_features
                                and canonical_feature not in self._only_features
                            ):
                                logging.warning(
                                    f"Feature {feature} not supported by quantms rescoring or not in only_features"
                                )
                            else:
                                oms_psm.setMetaValue(
                                    canonical_feature, OpenMSHelper.get_str_metavalue_round(value)
                                )
                                features.append(canonical_feature)
                        else:
                            logging.warning(
                                f"Feature {feature} not supported by quantms rescoring or not in only_features"
                            )
                hits.append(oms_psm)
            oms_peptide.setHits(hits)
            oms_peptides.append(oms_peptide)
        features = set(features)
        if features:
            logging.info(f"Features added to the peptides: {features} - adding them as meta value")
            search_parameters = self._idxml_reader.oms_proteins[0].getSearchParameters()
            try:
                features_existing = search_parameters.getMetaValue("extra_features")
            except Exception:
                logging.info("No extra features found in the search parameters")
                features_existing = None

            if features_existing is None:
                features_existing = ""
            extra_features = (features_existing + "," if features_existing else "") + ",".join(
                features
            )
            search_parameters.setMetaValue("extra_features", extra_features)
            self._idxml_reader.oms_proteins[0].setSearchParameters(search_parameters)

        self._idxml_reader.oms_peptides = oms_peptides

    def _get_top_batch_psms(self, psm_list: PSMList) -> PSMList:
        """
        Retrieve the top batch of PSMs for calibration based on their scores.

        This method sorts the provided PSMList by PSM score, taking into account
        whether a lower score is considered better. It then selects a subset of
        the sorted list to be used as a calibration set, determined by the
        `calibration_set_size` attribute.

        Parameters
        ----------
        psm_list : PSMList
            A list of peptide-spectrum matches to be sorted and filtered.

        Returns
        -------
        PSMList
            A subset of the input PSMList representing the top batch for calibration.
        """
        logging.info("Getting top calibration set from PSMs.")

        # Sort ms2pip results by PSM score and lower score is better
        ms2pip_results_copy = psm_list.psm_list.copy()

        ms2pip_results_copy = [result for result in ms2pip_results_copy if not result.is_decoy]
        ms2pip_results_copy.sort(key=lambda x: x.score, reverse=not self._lower_score_is_better)

        # Get a calibration set, the % of psms to be used for calibrarion is defined by calibration_set_size
        calibration_set = PSMList(
            psm_list=ms2pip_results_copy[
                : int(len(ms2pip_results_copy) * 0.6)  # Select the 60 of the PSMS.
            ]
        )

        return PSMList(psm_list=calibration_set)

    def _get_highest_fragmentation(self) -> Union[str, None]:
        """
        Determine the highest fragmentation method used in the dataset.

        This method retrieves the fragmentation statistics from the idXML reader
        and identifies the most frequently used fragmentation method, prioritizing
        "HCD" over "CID". If no statistics or methods are available, it logs a warning
        and returns None.

        Returns
        -------
        Union[str, None]
            The highest fragmentation method ("HCD" or "CID") or None if unavailable.
        """
        stats = self._idxml_reader.stats
        if stats is None or stats.ms_level_dissociation_method is None:
            logging.warning("No stats found or no ms_level_dissociation_methods found")
            return None

        first_key = max(
            stats.ms_level_dissociation_method, key=stats.ms_level_dissociation_method.get
        )
        if first_key[1] == "HCD":
            return "HCD"
        elif first_key[1] == "CID":
            return "CID"
        return None

    def _get_mae_from_psm_list(self, retrained_psms: PSMList) -> float:
        """
        Calculate the Mean Absolute Error (MAE) from a list of PSMs.

        This method calculates the Mean Absolute Error (MAE) from a list of PSMs
        by comparing the predicted retention time to the observed retention time
        for each PSM. The MAE is calculated as the average of the absolute differences
        between the predicted and observed retention times.

        Parameters
        ----------
        retrained_psms : PSMList
            A list of PSMs with predicted and observed retention times.

        Returns
        -------
        float
            The Mean Absolute Error (MAE) calculated from the PSMs.
        """
        best_scored_psms = self._get_top_batch_psms(retrained_psms)
        mae = 0
        for psm in best_scored_psms:
            mae += abs(psm.rescoring_features["rt_diff"])
        return mae / len(best_scored_psms)
