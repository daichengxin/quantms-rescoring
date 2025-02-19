from dataclasses import dataclass
import logging
from collections import defaultdict
from pathlib import Path
from typing import Union, Iterable, List, Optional, Dict, Tuple
from warnings import filterwarnings

from psm_utils import PSM, PSMList
import pyopenms as oms
from psm_utils.io.idxml import IdXMLReader

from quantmsrescore.openms import OpenMSHelper

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# Suppress OpenMS warning about data path
filterwarnings(
    "ignore",
    message="OPENMS_DATA_PATH environment variable already exists",
    category=UserWarning,
    module="pyopenms",
)


@dataclass
class ScoreStats:
    """Statistics about score occurrence in peptide hits."""

    total_hits: int
    missing_count: int

    @property
    def missing_percentage(self) -> float:
        """Calculate percentage of hits missing this score."""
        return (self.missing_count / self.total_hits * 100) if self.total_hits else 0


class IdXMLRescoringReader:
    """
    Reader class for processing and rescoring idXML files containing peptide identifications.

    This class handles reading and parsing idXML files, managing PSMs (Peptide-Spectrum Matches),
    and provides functionality for spectrum validation and scoring analysis.
    """

    def __init__(self, filename: Union[Path, str]) -> None:
        """
        Initialize the IdXMLRescoringReader.

        Args:
            filename: Path to the idXML file to process
        """

        self.filename = Path(filename)
        self.high_score_better: Optional[bool] = None
        self.skip_invalid_psm = 0
        self.new_peptide_ids: List[oms.PeptideIdentification] = []

        # Internal state
        self._mzml_path = None
        self._psms: PSMList = None
        self._spec_lookup: Optional[oms.SpectrumLookup] = None
        self._exp: Optional[oms.MSExperiment] = None

        # Parse the input file
        self.oms_proteins, self.oms_peptides = self._parse_idxml()

    def __iter__(self) -> Iterable[PSM]:
        """
        Iterate over valid PSMs in the idXML file.

        Yields:
            PSM objects for each valid peptide hit
        """
        score_stats = self._analyze_score_coverage()
        self._log_score_coverage(score_stats)

        for peptide_id in self.oms_peptides:
            valid_hits = []
            for peptide_hit in peptide_id.getHits():
                psm = self._parse_psm(self.oms_proteins, peptide_id, peptide_hit)
                if psm is not None:
                    valid_hits.append(peptide_hit)
                    yield psm
                else:
                    self.skip_invalid_psm += 1

            if valid_hits:
                peptide_id.setHits(valid_hits)
                self.new_peptide_ids.append(peptide_id)

    def _analyze_score_coverage(self) -> Dict[str, ScoreStats]:
        """
        Analyze the coverage of scores across all peptide hits.

        Returns:
            Dictionary mapping score names to their statistics
        """
        scores_stats = defaultdict(lambda: ScoreStats(0, 0))
        total_hits = 0

        for peptide_id in self.oms_peptides:
            hits = peptide_id.getHits()
            total_hits += len(hits)
            for hit in hits:
                meta_values = []
                hit.getKeys(meta_values)
                for score in meta_values:
                    scores_stats[score].total_hits += 1

        # Calculate missing counts
        for score, stats in scores_stats.items():
            stats.missing_count = total_hits - stats.total_hits

        return scores_stats

    def _log_score_coverage(self, score_stats: Dict[str, ScoreStats]) -> None:
        """
        Log warnings for scores with incomplete coverage.

        Args:
            score_stats: Dictionary of score statistics to analyze
        """
        for score, stats in score_stats.items():
            if stats.missing_count > 0:
                logging.warning(
                    f"Score {score} is missing in {stats.missing_count} PSMs "
                    f"({stats.missing_percentage:.1f}% of total)"
                )
                if stats.missing_percentage > 10:
                    logging.error(f"Score {score} is missing in more than 10% of PSMs")

    @staticmethod
    def _parse_psm(
        protein_ids: Union[oms.ProteinIdentification, List[oms.ProteinIdentification]],
        peptide_id: oms.PeptideIdentification,
        peptide_hit: oms.PeptideHit,
        is_decoy: bool = False,
    ) -> Optional[PSM]:
        """
        Parse OpenMS peptide hit data into a PSM object.

        Args:
            protein_ids: Protein identification data
            peptide_id: Peptide identification data
            peptide_hit: Individual peptide hit to parse
            is_decoy: Whether this hit is a decoy match

        Returns:
            Parsed PSM object, or None if parsing fails
        """
        try:
            peptidoform = IdXMLReader._parse_peptidoform(
                peptide_hit.getSequence().toString(), peptide_hit.getCharge()
            )

            # Create unique identifier for provenance tracking
            spectrum_ref = peptide_id.getMetaValue("spectrum_reference")
            key = f"{peptidoform}/{peptide_id.getRT()}/{spectrum_ref}"
            value = (
                f"{peptide_hit.getSequence().toString()}/"
                f"{peptide_hit.getCharge()}/{peptide_id.getRT()}/"
                f"{spectrum_ref}"
            )

            return PSM(
                peptidoform=peptidoform,
                spectrum_id=spectrum_ref,
                run=IdXMLReader._get_run(protein_ids, peptide_id),
                is_decoy=is_decoy,
                score=peptide_hit.getScore(),
                precursor_mz=peptide_id.getMZ(),
                retention_time=peptide_id.getRT(),
                rank=peptide_hit.getRank() + 1,  # Convert 0-based to 1-based
                source="idXML",
                provenance_data={key: value},
            )
        except Exception as e:
            logging.error(f"Failed to parse PSM: {e}")
            return None

    def _parse_idxml(
        self,
    ) -> Tuple[List[oms.ProteinIdentification], List[oms.PeptideIdentification]]:
        """
        Parse the idXML file to extract protein and peptide identifications.

        Returns:
            Tuple of (protein identifications, peptide identifications)
        """
        idxml_file = oms.IdXMLFile()
        proteins, peptides = [], []
        idxml_file.load(str(self.filename), proteins, peptides)
        return proteins, peptides

    def read_file(self) -> PSMList:
        """
        Read and parse all PSMs from the idXML file.

        Returns:
            List of parsed PSM objects
        """

        psm_list = []  # List to store parsed PSMs
        for peptide_id in self.oms_peptides:
            if self.high_score_better is None:
                self.high_score_better = peptide_id.isHigherScoreBetter()
            elif self.high_score_better != peptide_id.isHigherScoreBetter():
                logging.warning("Inconsistent score direction found in idXML file")

            for psm_hit in peptide_id.getHits():
                psm = self._parse_psm(
                    protein_ids=self.oms_proteins,
                    peptide_id=peptide_id,
                    peptide_hit=psm_hit,
                    is_decoy=OpenMSHelper.is_decoy_peptide_hit(psm_hit),
                )
                if psm is not None:
                    psm_list.append(psm)

        self._psms = PSMList(psm_list=psm_list)
        logging.info(f"Loaded {len(self._psms)} PSMs from {self.filename}")
        return self._psms

    def build_spectrum_lookup(self, mzml_file: Union[str, Path]) -> None:
        """
        Build a spectrum lookup index from an mzML file.

        Args:
            mzml_file: Path to the mzML file to index
        """
        self._mzml_path = str(mzml_file) if isinstance(mzml_file, Path) else mzml_file
        self._exp, self._spec_lookup = OpenMSHelper.get_spectrum_lookup_indexer(self._mzml_path)
        logging.info(f"Built SpectrumLookup from {self._mzml_path}")

    def validate_psm_spectrum_references(self) -> int:
        """
        Validate that all PSM spectrum references exist in the loaded spectra.

        Returns:
            Number of PSMs with missing spectra

        Raises:
            ValueError: If spectrum lookup or PSMs are not initialized
        """
        if self._spec_lookup is None or self._exp is None or not self._psms:
            raise ValueError("Spectrum lookup or PSMs not initialized")

        missing_count = 0
        ms_level_collections = defaultdict(int)
        for psm in self._psms:
            spectrum = OpenMSHelper.get_spectrum_for_psm(psm, self._exp, self._spec_lookup)
            if spectrum is None:
                logging.error(f"Spectrum not found for PSM {psm}")
                missing_count += 1
            num_peaks = len(spectrum.get_peaks()[0])
            if num_peaks == 0:
                logging.warning(f"Empty spectrum found for PSM {psm}")
                missing_count += 1
            ms_level = spectrum.getMSLevel()
            if ms_level == 3:
                logging.info(f"MS level 3 spectrum found for PSM {psm}")
            if ms_level > 0:
                ms_level_collections[ms_level] += 1
        logging.info(f"MS level distribution: {ms_level_collections}")

        if missing_count:
            logging.error(f"Found {missing_count} PSMs with missing spectra")

        return missing_count

    def get_spectrum_path(self):
        return self._mzml_path

    def get_psms(self):
        return self._psms

    def set_psms(self, psm_list):
        self._psms = psm_list
