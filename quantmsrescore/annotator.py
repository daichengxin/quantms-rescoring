import logging
from pathlib import Path
from typing import Union

from quantmsrescore.deeplc import DeepLCAnnotator
from quantmsrescore.idxmlreader import IdXMLRescoringReader
from quantmsrescore.ms2pip import MS2PIPAnnotator
from quantmsrescore.openms import OpenMSHelper

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


class Annotator:
    def __init__(
        self,
        feature_generators: str,
        ms2pip_model: str = "HCD2021",
        ms2pip_model_path: str = "models",
        ms2_tolerance: float = 0.05,
        calibration_set_size: float = 0.2,
        deeplc_retrain: bool = False,
        processes: int = 2,
        id_decoy_pattern: str = "^DECOY_",
        lower_score_is_better: bool = True,
        log_level: str = "INFO",
        spectrum_id_pattern: str = "(.*)",  # default for openms idXML
        psm_id_pattern: str = "(.*)",  # default for openms idXML
        remove_missing_spectra: bool = True,
        ms2_only: bool = True,
    ):
        """
        Initializes the Annotator class with configuration parameters for feature generation
        and rescoring of peptide-spectrum matches (PSMs).

        Parameters
        ----------
        feature_generators (str): Comma-separated list of feature generators to use (e.g., "deeplc,ms2pip").
        ms2pip_model (str): The MS2PIP model to use for annotation. Default is "HCD2021".
        ms2pip_model_path (str): Path to the directory containing MS2PIP models. Default is "models".
        ms2_tolerance (float): Tolerance for MS2PIP annotation. Default is 0.05.
        calibration_set_size (float): Fraction of data used for calibration. Default is 0.2.
        deeplc_retrain (bool): Whether to retrain DeepLC models. Default is False.
        processes (int): Number of processes to use for parallel computation. Default is 2.
        id_decoy_pattern (str): Regex pattern to identify decoy IDs. Default is "^DECOY_".
        lower_score_is_better (bool): Whether a lower score indicates a better match. Default is True.
        log_level (str): Logging level for the annotator. Default is "INFO".
        spectrum_id_pattern (str): Regex pattern for spectrum IDs. Default is "(.*)".
        psm_id_pattern (str): Regex pattern for PSM IDs. Default is "(.*)".

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
        self._deeplc_retrain = deeplc_retrain
        self._remove_missing_spectra = remove_missing_spectra
        self._ms2_only = ms2_only

    def build_idxml_data(self, idxml_file: Union[str, Path], spectrum_path: Union[str, Path]):
        """
        Build and load data from an idXML file for annotation.

        This method initializes the IdXMLRescoringReader with the specified
        idXML and mzML files, and counts the number of decoy and target PSMs.

        Parameters
        ----------
        idxml_file : Union[str, Path]
            The path to the idXML file to be processed.
        spectrum_path : Union[str, Path]
            The path to the corresponding mzML file.

        Logs
        ----
        Logs the number of PSMs, decoys, and targets loaded from the idXML file.
        """

        logging.info("Running the Annotator on file: %s", idxml_file)

        openms_helper = OpenMSHelper()

        # Load the idXML file and the corresponding mzML file
        self._idxml_reader = IdXMLRescoringReader(
            idexml_filename=idxml_file,
            mzml_file=spectrum_path,
            only_ms2=self._ms2_only,
            remove_missing_spectrum=self._remove_missing_spectra,
        )
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

        logging.debug(f"Running Annotations with following configurations: {self.__dict__}")

        if self._ms2pip:
            logging.info("Running MS2PIP on the PSMs")

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

            psm_list = self._idxml_reader.psms
            ms2pip_generator.add_features(psm_list)
            self._idxml_reader.psms = psm_list

            logging.info("MS2PIP Annotations added to the PSMs")

        if self._deepLC:
            logging.info("Running deepLC on the PSMs")

            deeplc_annotator = DeepLCAnnotator(
                self._lower_score_is_better,
                calibration_set_size=self._calibration_set_size,
                processes=self._processes,
            )
            psm_list = self._idxml_reader.psms
            deeplc_annotator.add_features(psm_list)
            self._idxml_reader.psms = psm_list

        if self._ms2pip or self._deepLC:
            self._convert_features_psms_to_oms_peptides()

        logging.info("Annotations added to the PSMs, starting to modified OMS peptides")

    def write_idxml_file(self, filename: Union[str, Path]):
        OpenMSHelper.write_idxml_file(
            filename=filename,
            protein_ids=self._idxml_reader.openms_proteins,
            peptide_ids=self._idxml_reader.openms_peptides,
        )
        logging.info("Annotated idXML file written to %s", filename)

    def _convert_features_psms_to_oms_peptides(self):
        """
        This method converts features from PSMs to OMS peptides features. This is required as OpenMS
        uses OMS peptides for storing features.
        """
        psm_dict = {next(iter(psm.provenance_data)): psm for psm in self._idxml_reader.psms}

        oms_peptides = []

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
                        # Round to 4 decimal str for OpenMS
                        value_str = "{:.4f}".format(value)
                        oms_psm.setMetaValue(feature, value_str)
                hits.append(oms_psm)
            oms_peptide.setHits(hits)
            oms_peptides.append(oms_peptide)
        self._idxml_reader.oms_peptides = oms_peptides
