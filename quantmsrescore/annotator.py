import copy
from pathlib import Path
from typing import Optional, Set, Union

from psm_utils import PSMList

from quantmsrescore.deeplc import DeepLCAnnotator
from quantmsrescore.exceptions import Ms2pipIncorrectModelException
from quantmsrescore.idxmlreader import IdXMLRescoringReader
from quantmsrescore.logging_config import get_logger
from quantmsrescore.ms2pip import MS2PIPAnnotator
from quantmsrescore.openms import OpenMSHelper

# Get logger for this module
logger = get_logger(__name__)


class FeatureAnnotator:
    """
    Annotator for peptide-spectrum matches (PSMs) using MS2PIP and DeepLC models.

    This class handles the annotation of PSMs with additional features generated
    from MS2PIP and DeepLC models to improve rescoring.
    """

    def __init__(
        self,
        feature_generators: str,
        only_features: Optional[str] = None,
        ms2pip_model: str = "HCD2021",
        force_model: bool = False,
        ms2pip_model_path: str = "models",
        ms2_tolerance: float = 0.05,
        calibration_set_size: float = 0.2,
        valid_correlations_size: float = 0.7,
        skip_deeplc_retrain: bool = False,
        processes: int = 2,
        log_level: str = "INFO",
        spectrum_id_pattern: str = "(.*)",  # default for openms idXML
        psm_id_pattern: str = "(.*)",  # default for openms idXML
        remove_missing_spectra: bool = True,
        ms2_only: bool = True
    ):
        """
        Initialize the Annotator with configuration parameters.

        Parameters
        ----------
        feature_generators : str
            Comma-separated list of feature generators (e.g., "ms2pip,deeplc").
        only_features : str, optional
            Comma-separated list of features to include in annotation.
        ms2pip_model : str, optional
            MS2PIP model name (default: "HCD2021").
        ms2pip_model_path : str, optional
            Path to MS2PIP model directory (default: "models").
        ms2_tolerance : float, optional
            MS2 tolerance for feature generation (default: 0.05).
        calibration_set_size : float, optional
            Percentage of PSMs to use for calibration (default: 0.2).
        skip_deeplc_retrain : bool, optional
            Skip retraining the deepLC model (default: False).
        processes : int, optional
            Number of parallel processes (default: 2).
        log_level : str, optional
            Logging level (default: "INFO").
        spectrum_id_pattern : str, optional
            Pattern for identifying spectrum IDs (default: "(.*)").
        psm_id_pattern : str, optional
            Pattern for identifying PSM IDs (default: "(.*)").
        remove_missing_spectra : bool, optional
            Remove PSMs with missing spectra (default: True).
        ms2_only : bool, optional
            Process only MS2-level PSMs (default: True).
        find_best_ms2pip_model : bool, optional
            Find best MS2PIP model for the dataset (default: False).

        Raises
        ------
        ValueError
            If no feature generators are provided or if neither ms2pip nor deeplc is specified.
        """
        # Set up logging
        from quantmsrescore.logging_config import configure_logging

        configure_logging(log_level)

        # Validate inputs
        if not feature_generators:
            raise ValueError("feature_generators must be provided.")

        feature_annotators = feature_generators.split(",")
        if not any(annotator in feature_annotators for annotator in ["deeplc", "ms2pip"]):
            raise ValueError("At least one of deeplc or ms2pip must be provided.")

        # Initialize state
        self._idxml_reader = None
        self._deepLC = "deeplc" in feature_annotators
        self._ms2pip = "ms2pip" in feature_annotators

        # Parse and validate features
        self._only_features = []
        if only_features:
            self._only_features = OpenMSHelper.validate_features(only_features.split(","))

        # Store configuration
        self._ms2pip_model = ms2pip_model
        self._ms2pip_model_path = ms2pip_model_path
        self._ms2_tolerance = ms2_tolerance
        self._calibration_set_size = calibration_set_size
        self._valid_correlations_size = valid_correlations_size
        self._processes = processes
        self._higher_score_better = None
        self._spectrum_id_pattern = spectrum_id_pattern
        self._psm_id_pattern = psm_id_pattern
        self._skip_deeplc_retrain = skip_deeplc_retrain
        self._remove_missing_spectra = remove_missing_spectra
        self._ms2_only = ms2_only
        self._force_model = force_model

    def build_idxml_data(
        self, idxml_file: Union[str, Path], spectrum_path: Union[str, Path]
    ) -> None:
        """
        Load data from idXML and mzML files.

        Parameters
        ----------
        idxml_file : Union[str, Path]
            Path to the idXML file containing PSM data.
        spectrum_path : Union[str, Path]
            Path to the corresponding mzML file with spectral data.

        Raises
        ------
        Exception
            If loading the files fails.
        """
        logger.info(f"Loading data from: {idxml_file}")

        try:
            # Convert paths to Path objects for consistency
            idxml_path = Path(idxml_file)
            spectrum_path = Path(spectrum_path)

            # Load the idXML file and corresponding mzML file
            self._idxml_reader = IdXMLRescoringReader(
                idxml_filename=idxml_path,
                mzml_file=spectrum_path,
                only_ms2=self._ms2_only,
                remove_missing_spectrum=self._remove_missing_spectra,
            )
            self._higher_score_better = self._idxml_reader.high_score_better

            # Log statistics about loaded data
            psm_list = self._idxml_reader.psms
            openms_helper = OpenMSHelper()
            decoys, targets = openms_helper.count_decoys_targets(self._idxml_reader.oms_peptides)

            logger.info(
                f"Loaded {len(psm_list)} PSMs from {idxml_path.name}: {decoys} decoys and {targets} targets"
            )

        except Exception as e:
            logger.error(f"Failed to load input files: {str(e)}")
            raise

    def annotate(self) -> None:
        """
        Annotate PSMs with MS2PIP and/or DeepLC features.

        This method runs the selected feature generators to add annotations
        to the loaded PSMs.

        Raises
        ------
        ValueError
            If no idXML data is loaded.
        """
        if not self._idxml_reader:
            raise ValueError("No idXML data loaded. Call build_idxml_data() first.")

        logger.debug(f"Running annotations with configuration: {self.__dict__}")

        # Run MS2PIP annotation if enabled
        if self._ms2pip:
            self._run_ms2pip_annotation()

        # Run DeepLC annotation if enabled
        if self._deepLC:
            self._run_deeplc_annotation()

        # Convert features to OpenMS format if any annotations were added
        if self._ms2pip or self._deepLC:
            self._convert_features_psms_to_oms_peptides()

        logger.info("Annotation complete")

    def write_idxml_file(self, filename: Union[str, Path]) -> None:
        """
        Write annotated data to idXML file.

        Parameters
        ----------
        filename : Union[str, Path]
            Path where the annotated idXML file will be written.

        Raises
        ------
        Exception
            If writing the file fails.
        """
        try:
            out_path = Path(filename)
            OpenMSHelper.write_idxml_file(
                filename=out_path,
                protein_ids=self._idxml_reader.openms_proteins,
                peptide_ids=self._idxml_reader.openms_peptides,
            )
            logger.info(f"Annotated idXML file written to {out_path}")
        except Exception as e:
            logger.error(f"Failed to write annotated idXML file: {str(e)}")
            raise

    def _run_ms2pip_annotation(self) -> None:
        """Run MS2PIP annotation on the loaded PSMs."""
        logger.info("Running MS2PIP annotation")

        # Initialize MS2PIP annotator
        try:
            ms2pip_generator = self._create_ms2pip_annotator()
        except Exception as e:
            logger.error(f"Failed to initialize MS2PIP: {e}")
            raise

        # Apply MS2PIP annotation
        psm_list = self._idxml_reader.psms
        try:
            ms2pip_generator.add_features(psm_list)
            self._idxml_reader.psms = psm_list
            logger.info("MS2PIP annotations added to PSMs")
        except Ms2pipIncorrectModelException:
            if not self._force_model:
                self._find_and_apply_best_ms2pip_model(psm_list)
            else:
                logger.error("MS2PIP model not suitable for this data. Skip MS2PIP annotations.")

        except Exception as e:
            logger.error(f"Failed to add MS2PIP features: {e}")

    def _create_ms2pip_annotator(
        self, model: Optional[str] = None, tolerance: Optional[float] = None
    ) -> MS2PIPAnnotator:
        """
        Create an MS2PIP annotator with the specified or default model.

        Parameters
        ----------
        model : str, optional
            MS2PIP model name to use, defaults to self._ms2pip_model if None.

        Returns
        -------
        MS2PIPAnnotator
            Configured MS2PIP annotator.
        """
        return MS2PIPAnnotator(
            ms2_tolerance=tolerance or self._ms2_tolerance,
            model=model or self._ms2pip_model,
            spectrum_path=self._idxml_reader.spectrum_path,
            spectrum_id_pattern=self._spectrum_id_pattern,
            model_dir=self._ms2pip_model_path,
            calibration_set_size=self._calibration_set_size,
            valid_correlations_size=self._valid_correlations_size,
            correlation_threshold=0.7,  # Consider making this configurable
            higher_score_better=self._higher_score_better,
            processes=self._processes,
            force_model=self._force_model
        )

    def _find_and_apply_best_ms2pip_model(self, psm_list: PSMList) -> None:
        """
        Find and apply the best MS2PIP model for the dataset.

        Parameters
        ----------
        psm_list : PSMList
            List of PSMs to annotate.
        """
        logger.info("Finding best MS2PIP model for the dataset")

        # Get top scoring PSMs for model selection
        batch_psms = self._get_top_batch_psms(psm_list)

        # Create annotator with default model to use for finding best model
        ms2pip_generator = self._create_ms2pip_annotator()

        # Find best model based on fragmentation type
        fragmentation = self._get_highest_fragmentation()
        model, corr = ms2pip_generator._find_best_ms2pip_model(
            batch_psms=batch_psms,
            knwon_fragmentation=fragmentation,
        )

        if model:
            logger.info(f"Best model found: {model} with average correlation {corr}")

            # Create new annotator with best model
            ms2pip_generator = self._create_ms2pip_annotator(model=model, tolerance=self._ms2_tolerance)

            # Apply annotation with best model
            ms2pip_generator.add_features(psm_list)
            self._idxml_reader.psms = psm_list
            logger.info("MS2PIP annotations added using best model")
        else:
            logger.error("No suitable MS2PIP model found for this dataset. Skip MS2PIP annotations.")

    def _run_deeplc_annotation(self) -> None:
        """Run DeepLC annotation on the loaded PSMs."""
        logger.info("Running DeepLC annotation")

        try:
            if self._skip_deeplc_retrain:
                # Simple case - use pre-trained model
                deeplc_annotator = self._create_deeplc_annotator(retrain=False)
            else:
                # Compare retrained vs pretrained performance
                deeplc_annotator = self._determine_optimal_deeplc_model()

            # Apply annotation
            psm_list = self._idxml_reader.psms
            deeplc_annotator.add_features(psm_list)
            self._idxml_reader.psms = psm_list
            logger.info("DeepLC annotations added to PSMs")

        except Exception as e:
            logger.error(f"Failed to apply DeepLC annotation: {e}")
            raise

    def _create_deeplc_annotator(
        self, retrain: bool = False, calibration_set_size: float = None
    ) -> DeepLCAnnotator:
        """
        Create a DeepLC annotator with specified configuration.

        Parameters
        ----------
        retrain : bool
            Whether to retrain the DeepLC model.

        Returns
        -------
        DeepLCAnnotator
            Configured DeepLC annotator.
        """
        kwargs = {"deeplc_retrain": retrain}

        if calibration_set_size is None:
            calibration_set_size = self._calibration_set_size

        return DeepLCAnnotator(
            not self._higher_score_better,
            calibration_set_size=calibration_set_size,
            processes=self._processes,
            **kwargs,
        )

    def _determine_optimal_deeplc_model(self) -> DeepLCAnnotator:
        """
        Determine the optimal DeepLC model by comparing retrained vs. pretrained performance.

        This function evaluates both a retrained model and a pretrained model on the same dataset,
        calculates the Mean Absolute Error (MAE) for each, and selects the model with lower error.

        Returns
        -------
        DeepLCAnnotator
            The DeepLC annotator with the lowest MAE (best performance).
        """
        # Get base PSMs for comparison
        base_psms = self._idxml_reader.psms.psm_list

        # Evaluate retrained model
        retrained_psms = PSMList(psm_list=copy.deepcopy(base_psms))
        retrained_model = self._create_deeplc_annotator(retrain=True, calibration_set_size=0.6)
        retrained_model.add_features(retrained_psms)
        mae_retrained = self._get_mae_from_psm_list(retrained_psms)

        # Evaluate pretrained model
        pretrained_psms = PSMList(psm_list=copy.deepcopy(base_psms))
        pretrained_model = self._create_deeplc_annotator(retrain=False, calibration_set_size=0.6)
        pretrained_model.add_features(pretrained_psms)
        mae_pretrained = self._get_mae_from_psm_list(pretrained_psms)

        # Select model with lower MAE
        if mae_retrained < mae_pretrained:
            logger.info(
                f"Retrained DeepLC model has lower MAE ({mae_retrained:.4f} vs {mae_pretrained:.4f}), using it: {retrained_model.selected_model}"
            )
            return retrained_model
        else:
            logger.info(
                f"Pretrained DeepLC model has lower/equal MAE ({mae_pretrained:.4f} vs {mae_retrained:.4f}), using it: {pretrained_model.selected_model}"
            )
            return pretrained_model

    def _convert_features_psms_to_oms_peptides(self) -> None:
        """
        Transfer features from PSM objects to OpenMS peptide objects.
        """
        # Create lookup dictionary for PSMs
        psm_dict = {next(iter(psm.provenance_data)): psm for psm in self._idxml_reader.psms}

        oms_peptides = []
        added_features: Set[str] = set()

        # Process each peptide
        for oms_peptide in self._idxml_reader.oms_peptides:
            hits = []

            # Process each hit within the peptide
            for oms_psm in oms_peptide.getHits():
                psm_hash = OpenMSHelper.get_psm_hash_unique_id(
                    peptide_hit=oms_peptide, psm_hit=oms_psm
                )

                psm = psm_dict.get(psm_hash)
                if psm is None:
                    logger.warning(f"PSM not found for peptide {oms_peptide.getMetaValue('id')}")
                else:
                    # Add features to the OpenMS PSM
                    for feature, value in psm.rescoring_features.items():
                        canonical_feature = OpenMSHelper.get_canonical_feature(feature)

                        if canonical_feature is not None:
                            if (
                                self._only_features
                                and canonical_feature not in self._only_features
                            ):
                                continue

                            oms_psm.setMetaValue(
                                canonical_feature, OpenMSHelper.get_str_metavalue_round(value)
                            )
                            added_features.add(canonical_feature)
                        else:
                            logger.debug(f"Feature {feature} not supported by quantms rescoring")

                hits.append(oms_psm)

            oms_peptide.setHits(hits)
            oms_peptides.append(oms_peptide)

        # Update search parameters with added features
        self._update_search_parameters(added_features)

        # Update the peptides in the reader
        self._idxml_reader.oms_peptides = oms_peptides

    def _update_search_parameters(self, features: Set[str]) -> None:
        """
        Update search parameters with new features.

        Parameters
        ----------
        features : Set[str]
            Set of feature names to add to search parameters.
        """
        if not features:
            return

        logger.info(f"Adding features to search parameters: {', '.join(sorted(features))}")

        # Get search parameters
        search_parameters = self._idxml_reader.oms_proteins[0].getSearchParameters()

        # Get existing features
        try:
            features_existing = search_parameters.getMetaValue("extra_features")
            if features_existing:
                existing_set = set(features_existing.split(","))
            else:
                existing_set = set()
        except Exception:
            existing_set = set()

        # Combine existing and new features
        all_features = existing_set.union(features)

        # Update search parameters
        search_parameters.setMetaValue("extra_features", ",".join(sorted(all_features)))
        self._idxml_reader.oms_proteins[0].setSearchParameters(search_parameters)

    def _get_top_batch_psms(self, psm_list: PSMList) -> PSMList:
        """
        Get top-scoring non-decoy PSMs for calibration.

        Parameters
        ----------
        psm_list : PSMList
            List of PSMs to filter.

        Returns
        -------
        PSMList
            Filtered list containing top-scoring PSMs.
        """
        logger.info("Selecting top PSMs for calibration")

        # Filter non-decoy PSMs
        non_decoy_psms = [result for result in psm_list.psm_list if not result.is_decoy]

        if not non_decoy_psms:
            logger.warning("No non-decoy PSMs found for calibration")
            return PSMList(psm_list=[])

        # Sort by score
        non_decoy_psms.sort(key=lambda x: x.score, reverse=self._higher_score_better)

        # Select top 60% for calibration
        calibration_size = max(1, int(len(non_decoy_psms) * 0.6))
        calibration_psms = non_decoy_psms[:calibration_size]

        return PSMList(psm_list=calibration_psms)

    def _get_highest_fragmentation(self) -> Optional[str]:
        """
        Determine the predominant fragmentation method in the dataset.

        Returns
        -------
        Optional[str]
            "HCD", "CID", or None if not determined.
        """
        stats = self._idxml_reader.stats
        if not stats or not stats.ms_level_dissociation_method:
            logger.warning("No fragmentation method statistics available")
            return None

        # Find the most common fragmentation method
        most_common = max(
            stats.ms_level_dissociation_method, key=stats.ms_level_dissociation_method.get
        )

        # Return "HCD" or "CID" if applicable
        if most_common[1] in ["HCD", "CID"]:
            return most_common[1]

        return None

    def _get_mae_from_psm_list(self, psm_list: PSMList) -> float:
        """
        Calculate Mean Absolute Error of retention time prediction.

        Parameters
        ----------
        psm_list : PSMList
            List of PSMs with retention time predictions.

        Returns
        -------
        float
            Mean Absolute Error (MAE) value or infinity if calculation fails.
        """
        best_scored_psms = self._get_top_batch_psms(psm_list)

        if not best_scored_psms.psm_list:
            logger.warning("No PSMs available for MAE calculation")
            return float("inf")

        total_error = 0.0
        count = 0

        for psm in best_scored_psms.psm_list:
            if "rt_diff" in psm.rescoring_features:
                total_error += abs(psm.rescoring_features["rt_diff"])
                count += 1

        if count == 0:
            logger.warning("No valid retention time differences for MAE calculation")
            return float("inf")

        return total_error / count
