# Get logger for this module
from quantmsrescore.logging_config import get_logger

logger = get_logger(__name__)

from collections import defaultdict
from pathlib import Path
from typing import Union, List, Optional, Dict, Tuple, DefaultDict
from warnings import filterwarnings

filterwarnings(
    "ignore",
    message="OPENMS_DATA_PATH environment variable already exists",
    category=UserWarning,
    module="pyopenms",
)

import psm_utils
import pyopenms as oms
from psm_utils import PSM, PSMList
from pyopenms import IDFilter

from quantmsrescore.exceptions import MS3NotSupportedException
from quantmsrescore.openms import OpenMSHelper


class ScoreStats:
    """Statistics about score occurrence in peptide hits."""

    def __init__(self):
        self.total_hits: int = 0
        self.missing_count: int = 0

    @property
    def missing_percentage(self) -> float:
        """Calculate percentage of missing scores."""
        return (self.missing_count / self.total_hits * 100) if self.total_hits else 0


class SpectrumStats:
    """Statistics about spectrum analysis."""

    def __init__(self):
        self.missing_spectra: int = 0
        self.empty_spectra: int = 0
        self.ms_level_counts: DefaultDict[int, int] = defaultdict(int)
        self.ms_level_dissociation_method: Dict[Tuple[int, str], int] = {}
        self.predicted_ms_tolerance: Tuple[float, Optional[str]] = (0.0, None)
        self.reported_ms_tolerance: Tuple[float, Optional[str]] = (0.0, None)


class IdXMLReader:
    """
    A class to read and parse idXML files for protein and peptide identifications.

    Attributes
    ----------
    filename : Path
        The path to the idXML file.
    oms_proteins : List[oms.ProteinIdentification]
        List of protein identifications parsed from the idXML file.
    oms_peptides : List[oms.PeptideIdentification]
        List of peptide identifications parsed from the idXML file.
    """

    def __init__(self, idexml_filename: Union[Path, str]) -> None:
        """
        Initialize IdXMLReader with the specified idXML file.

        Parameters
        ----------
        idexml_filename : Union[Path, str]
            Path to the idXML file to be read and parsed.
        """
        self.filename = Path(idexml_filename)
        self.oms_proteins, self.oms_peptides = self._parse_idxml()

        # Private properties for spectrum lookup
        self._spec_lookup = None
        self._exp = None
        self._mzml_path = None
        self._stats = None  # IdXML stats

    def _parse_idxml(
        self,
    ) -> Tuple[List[oms.ProteinIdentification], List[oms.PeptideIdentification]]:
        """
        Parse the idXML file to extract protein and peptide identifications.

        Returns
        -------
        Tuple[List[oms.ProteinIdentification], List[oms.PeptideIdentification]]
            A tuple containing lists of protein and peptide identifications.
        """
        idxml_file = oms.IdXMLFile()
        proteins, peptides = [], []
        idxml_file.load(str(self.filename), proteins, peptides)
        return proteins, peptides

    @property
    def openms_proteins(self) -> List[oms.ProteinIdentification]:
        """Get the list of protein identifications."""
        return self.oms_proteins

    @property
    def openms_peptides(self) -> List[oms.PeptideIdentification]:
        """Get the list of peptide identifications."""
        return self.oms_peptides

    @property
    def stats(self) -> Optional[SpectrumStats]:
        """Get spectrum statistics."""
        return self._stats

    @property
    def spectrum_path(self) -> Optional[Union[str, Path]]:
        """Get the path to the mzML file."""
        return self._mzml_path

    def build_spectrum_lookup(
        self, mzml_file: Union[str, Path], check_unix_compatibility: bool = False
    ) -> None:
        """
        Build a SpectrumLookup indexer from an mzML file.

        Parameters
        ----------
        mzml_file : Union[str, Path]
            The path to the mzML file to be processed.
        check_unix_compatibility : bool, optional
            Flag to check for Unix compatibility in the mzML file, by default, False.
        """
        self._mzml_path = str(mzml_file) if isinstance(mzml_file, Path) else mzml_file
        if check_unix_compatibility:
            OpenMSHelper.check_unix_compatibility(self._mzml_path)
        self._exp, self._spec_lookup = OpenMSHelper.get_spectrum_lookup_indexer(self._mzml_path)
        logger.info(f"Built SpectrumLookup from {self._mzml_path}")


class IdXMLRescoringReader(IdXMLReader):
    """
    Reader class for processing and rescoring idXML files containing peptide identifications.

    This class handles reading and parsing idXML files, managing PSMs (Peptide-Spectrum Matches),
    and provides functionality for spectrum validation and scoring analysis.

    Attributes
    ----------
    filename : Path
        Path to the idXML file.
    high_score_better : Optional[bool]
        Indicates if higher scores are better.
    """

    def __init__(
        self,
        idexml_filename: Union[Path, str],
        mzml_file: Union[str, Path],
        only_ms2: bool = True,
        remove_missing_spectrum: bool = True,
    ) -> None:
        """
        Initialize the IdXMLRescoringReader with the specified files.

        Parameters
        ----------
        idexml_filename : Union[Path, str]
            Path to the idXML file to be read and parsed.
        mzml_file : Union[str, Path]
            Path to the mzML file for spectrum lookup.
        only_ms2 : bool, optional
            Flag to filter for MS2 spectra only, by default True.
        remove_missing_spectrum : bool, optional
            Flag to remove PSMs with missing spectra, by default True.
        """
        super().__init__(idexml_filename)
        self.high_score_better: Optional[bool] = None

        # Private attributes
        self._psms: Optional[PSMList] = None
        self.build_spectrum_lookup(mzml_file, check_unix_compatibility=True)
        self._validate_psm_spectrum_references(
            only_ms2=only_ms2, remove_missing_spectrum=remove_missing_spectrum
        )
        self._build_psm_index(only_ms2=only_ms2)

    @property
    def psms(self) -> Optional[PSMList]:
        """Get the list of PSMs."""
        return self._psms

    @psms.setter
    def psms(self, psm_list: PSMList) -> None:
        """Set the list of PSMs."""
        if not isinstance(psm_list, PSMList):
            raise TypeError("psm_list must be an instance of PSMList")
        self._psms = psm_list

    def analyze_score_coverage(self) -> Dict[str, ScoreStats]:
        """
        Analyze the coverage of scores across peptide hits.

        Returns
        -------
        Dict[str, ScoreStats]
            A dictionary mapping score names to their respective statistics.
        """
        scores_stats: Dict[str, ScoreStats] = defaultdict(ScoreStats)
        total_hits = sum(len(peptide_id.getHits()) for peptide_id in self.oms_peptides)

        for peptide_id in self.oms_peptides:
            for hit in peptide_id.getHits():
                meta_values = []
                hit.getKeys(meta_values)
                for score in meta_values:
                    scores_stats[score].total_hits += 1

        for stats in scores_stats.values():
            stats.missing_count = total_hits - stats.total_hits

        return scores_stats

    @staticmethod
    def log_score_coverage(score_stats: Dict[str, ScoreStats]) -> None:
        """
        Log information about score coverage.

        Parameters
        ----------
        score_stats : Dict[str, ScoreStats]
            Dictionary mapping score names to their statistics.
        """
        for score, stats in score_stats.items():
            if stats.missing_count > 0:
                percentage = stats.missing_percentage
                logger.warning(
                    f"Score {score} is missing in {stats.missing_count} PSMs "
                    f"({percentage:.1f}% of total)"
                )
                if percentage > 10:
                    logger.error(f"Score {score} is missing in more than 10% of PSMs")

    @staticmethod
    def _parse_psm(
        protein_ids: Union[oms.ProteinIdentification, List[oms.ProteinIdentification]],
        peptide_id: oms.PeptideIdentification,
        peptide_hit: oms.PeptideHit,
        is_decoy: bool = False,
    ) -> Optional[PSM]:
        """
        Parse a peptide-spectrum match (PSM) from given protein and peptide models.

        Parameters
        ----------
        protein_ids : Union[oms.ProteinIdentification, List[oms.ProteinIdentification]]
            Protein identification(s) associated with the PSM.
        peptide_id : oms.PeptideIdentification
            Peptide identification containing the peptide hit.
        peptide_hit : oms.PeptideHit
            Peptide hit to be parsed into a PSM.
        is_decoy : bool, optional
            Indicates if the PSM is a decoy, by default False.

        Returns
        -------
        Optional[PSM]
            A PSM object if parsing is successful, otherwise None.
        """
        try:
            peptidoform = psm_utils.io.idxml.IdXMLReader._parse_peptidoform(
                peptide_hit.getSequence().toString(), peptide_hit.getCharge()
            )

            spectrum_ref = peptide_id.getMetaValue("spectrum_reference")
            rt = peptide_id.getRT()

            # Create provenance tracking models
            provenance_key = OpenMSHelper.get_psm_hash_unique_id(
                peptide_hit=peptide_id, psm_hit=peptide_hit
            )

            return PSM(
                peptidoform=peptidoform,
                spectrum_id=spectrum_ref,
                run=psm_utils.io.idxml.IdXMLReader._get_run(protein_ids, peptide_id),
                is_decoy=is_decoy,
                score=peptide_hit.getScore(),
                precursor_mz=peptide_id.getMZ(),
                retention_time=rt,
                rank=peptide_hit.getRank() + 1,  # Ranks in idXML start at 0
                source="idXML",
                provenance_data={provenance_key: ""},  # We use only the key for provenance
            )
        except Exception as e:
            logger.error(f"Failed to parse PSM: {e}")
            return None

    def _build_psm_index(self, only_ms2: bool = True) -> PSMList:
        """
        Read and parse the idXML file to extract PSMs.

        Parameters
        ----------
        only_ms2 : bool, optional
            Flag to filter for MS2 spectra only, by default True.

        Returns
        -------
        PSMList
            A list of parsed PSM objects.
        """
        psm_list = []

        if only_ms2 and self._spec_lookup is None:
            logger.warning("Spectrum lookup not initialized, cannot filter for MS2 spectra")
            only_ms2 = False

        for peptide_id in self.oms_peptides:
            if self.high_score_better is None:
                self.high_score_better = peptide_id.isHigherScoreBetter()
            elif self.high_score_better != peptide_id.isHigherScoreBetter():
                logger.warning("Inconsistent score direction found in idXML file")

            for psm_hit in peptide_id.getHits():
                if (
                    only_ms2
                    and self._spec_lookup is not None
                    and OpenMSHelper.get_ms_level(peptide_id, self._spec_lookup, self._exp) != 2
                ):
                    continue
                psm = self._parse_psm(
                    protein_ids=self.oms_proteins,
                    peptide_id=peptide_id,
                    peptide_hit=psm_hit,
                    is_decoy=OpenMSHelper.is_decoy_peptide_hit(psm_hit),
                )
                if psm is not None:
                    psm_list.append(psm)

        self._psms = PSMList(psm_list=psm_list)
        logger.info(f"Loaded {len(self._psms)} PSMs from {self.filename}")
        return self._psms

    def _validate_psm_spectrum_references(
        self, remove_missing_spectrum: bool = True, only_ms2: bool = True
    ) -> SpectrumStats:
        """
        Validate spectrum references for peptide identifications and filter based on criteria.

        This method validates each peptide identification by checking if its referenced
        spectrum exists and has peaks. It also tracks MS level statistics and dissociation
        methods. Optionally, removes peptide identifications with missing/empty spectra or
        those that are not MS2 level.

        Parameters
        ----------
        remove_missing_spectrum : bool, optional
            If True, removes peptide identifications with missing or empty spectra,
            by default True.
        only_ms2 : bool, optional
            If True, removes peptide identifications that reference non-MS2 spectra,
            by default True.

        Returns
        -------
        SpectrumStats
            Statistics about spectrum validation including counts of missing spectra,
            empty spectra, MS level distribution, and dissociation methods.

        Raises
        ------
        ValueError
            If spectrum lookup or experiment are not initialized.
        MS3NotSupportedException
            If MS3 spectra are found while only_ms2 is True.

        Notes
        -----
        This method modifies the internal list of peptide identifications by filtering
        out entries that don't meet the specified criteria. It also updates protein
        identifications to remove entries that no longer have associated peptides.
        """

        if self._spec_lookup is None or self._exp is None:
            raise ValueError("Spectrum lookup or PSMs not initialized")

        self._stats = SpectrumStats()

        new_peptide_ids = []
        peptide_removed = 0

        for peptide_id in self.oms_peptides:
            spectrum = OpenMSHelper.get_spectrum_for_psm(peptide_id, self._exp, self._spec_lookup)
            spectrum_reference = OpenMSHelper.get_spectrum_reference(peptide_id)

            missing_spectrum, empty_spectrum = False, False
            ms_level = 2

            if spectrum is None:

                logger.error(
                    f"Spectrum not found for PeptideIdentification with {spectrum_reference}"
                )
                self._stats.missing_spectra += 1
                missing_spectrum = True
            else:
                peaks = spectrum.get_peaks()[0]
                if peaks is None or len(peaks) == 0:
                    logger.warning(f"Empty spectrum found for PSM {spectrum_reference}")
                    empty_spectrum = True
                    self._stats.empty_spectra += 1

                ms_level = spectrum.getMSLevel()
                self._stats.ms_level_counts[ms_level] += 1

                self._process_dissociation_methods(spectrum, ms_level)

                if ms_level != 2 and only_ms2:
                    logger.info(
                        f"MS level {ms_level} spectrum found for PSM {spectrum_reference}. "
                        "MS2pip models are not trained on MS3 spectra"
                    )

            if (remove_missing_spectrum and (missing_spectrum or empty_spectrum)) or (
                only_ms2 and ms_level != 2
            ):
                logger.debug(f"Removing PSM {spectrum_reference}")
                peptide_removed += 1
            else:
                new_peptide_ids.append(peptide_id)

        if peptide_removed > 0:
            logger.warning(
                f"Removed {peptide_removed} PSMs with missing or empty spectra or MS3 spectra"
            )
            self.oms_peptides = new_peptide_ids
            oms_filter = IDFilter()
            # We only want to have protein accessions with at least one peptide identification
            oms_filter.removeEmptyIdentifications(self.oms_peptides)
            oms_filter.removeUnreferencedProteins(self.oms_proteins, self.oms_peptides)

        ms_tolerance, ms_unit = OpenMSHelper.get_ms_tolerance(self.oms_proteins)
        self._stats.reported_ms_tolerance = (ms_tolerance, ms_unit)

        self._log_spectrum_statistics()

        if only_ms2 and self._stats.ms_level_counts.get(3, 0) > 0:
            ms2_dissociation_methods = self._stats.ms_level_dissociation_method.get((2, "HCD"), 0)
            logger.error(
                "MS3 spectra found in MS2-only mode, please filter your search for MS2 or dissociation method: {}".format(
                    ms2_dissociation_methods
                )
            )
            raise MS3NotSupportedException("MS3 spectra found in MS2-only mode")

        if self._stats.reported_ms_tolerance[1] == "ppm":
            self._stats.predicted_ms_tolerance = OpenMSHelper.get_predicted_ms_tolerance(
                exp=self._exp, ppm_tolerance=self._stats.reported_ms_tolerance[0]
            )

        return self._stats

    def _process_dissociation_methods(self, spectrum, ms_level):
        """Process dissociation methods from spectrum precursors."""
        oms_dissociation_matrix = OpenMSHelper.get_pyopenms_dissociation_matrix()
        for precursor in spectrum.getPrecursors():
            for method_index in precursor.getActivationMethods():
                if (oms_dissociation_matrix is not None) and (
                    0 <= method_index < len(oms_dissociation_matrix)
                ):
                    method = (
                        ms_level,
                        OpenMSHelper.get_dissociation_method(
                            method_index, oms_dissociation_matrix
                        ),
                    )
                    self._stats.ms_level_dissociation_method[method] = (
                        self._stats.ms_level_dissociation_method.get(method, 0) + 1
                    )
                else:
                    logger.warning(f"Unknown dissociation method index {method_index}")

    def _log_spectrum_statistics(self):
        """Log statistics about spectrum validation."""
        if self._stats.missing_spectra or self._stats.empty_spectra:
            logger.error(
                f"Found {self._stats.missing_spectra} PSMs with missing spectra and "
                f"{self._stats.empty_spectra} PSMs with empty spectra"
            )

        if len({k[1] for k in self._stats.ms_level_dissociation_method}) > 1:
            logger.error(
                "Found multiple dissociation methods in the same MS level. "
                "MS2pip models are not trained for multiple dissociation methods"
            )

        logger.info(f"MS level distribution: {dict(self._stats.ms_level_counts)}")
        logger.info(
            f"Dissociation Method Distribution: {self._stats.ms_level_dissociation_method}"
        )
