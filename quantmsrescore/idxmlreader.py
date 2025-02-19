from dataclasses import dataclass
import logging
from collections import defaultdict
from pathlib import Path
from typing import Union, Iterable, List, Optional, Dict, Tuple, DefaultDict
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


class ScoreStats:
    """Statistics about score occurrence in peptide hits."""
    total_hits: int = 0
    missing_count: int = 0

    @property
    def missing_percentage(self) -> float:
        """Calculate percentage of hits missing this score."""
        return (self.missing_count / self.total_hits * 100) if self.total_hits else 0


class SpectrumStats:
    """Statistics about spectrum analysis."""
    missing_spectra: int = 0
    empty_spectra: int = 0
    ms_level_counts: DefaultDict[int, int] = defaultdict(int)


class IdXMLRescoringReader:
    """
    Reader class for processing and rescoring idXML files containing peptide identifications.

    This class handles reading and parsing idXML files, managing PSMs (Peptide-Spectrum Matches),
    and provides functionality for spectrum validation and scoring analysis.

    Attributes:
        filename (Path): Path to the idXML file
        high_score_better (Optional[bool]): Indicates if higher scores are better
        skip_invalid_psm (int): Counter for skipped invalid PSMs
    """

    def __init__(self, filename: Union[Path, str]) -> None:
        """Initialize the IdXMLRescoringReader."""
        self.filename = Path(filename)
        self.high_score_better: Optional[bool] = None
        self.skip_invalid_psm: int = 0
        self.new_peptide_ids: List[oms.PeptideIdentification] = []

        # Private attributes
        self._mzml_path: Optional[str] = None
        self._psms: Optional[PSMList] = None
        self._spec_lookup: Optional[oms.SpectrumLookup] = None
        self._exp: Optional[oms.MSExperiment] = None

        # Parse input file
        self.oms_proteins, self.oms_peptides = self._parse_idxml()

    @property
    def spectrum_path(self) -> Union[str | Path]:
        """Get the path to the spectrum file."""
        return self._mzml_path

    @property
    def psms(self) -> Optional[PSMList]:
        """Get the current PSM list."""
        return self._psms

    @psms.setter
    def psms(self, psm_list: PSMList) -> None:
        """Set the PSM list."""
        if not isinstance(psm_list, PSMList):
            raise TypeError("psm_list must be an instance of PSMList")
        self._psms = psm_list

    @property
    def openms_proteins(self) -> List[oms.ProteinIdentification]:
        """Get the OpenMS protein identifications."""
        return self.oms_proteins

    @property
    def openms_peptides(self) -> List[oms.PeptideIdentification]:
        """Get the OpenMS peptide identifications."""
        return self.oms_peptides

    def __iter__(self) -> Iterable[PSM]:
        """Iterate over valid PSMs in the idXML file."""
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
        """Analyze the coverage of scores across all peptide hits."""
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

    def _log_score_coverage(self, score_stats: Dict[str, ScoreStats]) -> None:
        """Log warnings for scores with incomplete coverage."""
        for score, stats in score_stats.items():
            if stats.missing_count > 0:
                percentage = stats.missing_percentage
                logging.warning(
                    f"Score {score} is missing in {stats.missing_count} PSMs "
                    f"({percentage:.1f}% of total)"
                )
                if percentage > 10:
                    logging.error(f"Score {score} is missing in more than 10% of PSMs")

    @staticmethod
    def _parse_psm(
            protein_ids: Union[oms.ProteinIdentification, List[oms.ProteinIdentification]],
            peptide_id: oms.PeptideIdentification,
            peptide_hit: oms.PeptideHit,
            is_decoy: bool = False,
    ) -> Optional[PSM]:
        """Parse OpenMS peptide hit data into a PSM object."""
        try:
            peptidoform = IdXMLReader._parse_peptidoform(
                peptide_hit.getSequence().toString(),
                peptide_hit.getCharge()
            )

            spectrum_ref = peptide_id.getMetaValue("spectrum_reference")
            rt = peptide_id.getRT()

            # Create provenance tracking data
            provenance_key = f"{peptidoform}/{rt}/{spectrum_ref}"
            provenance_value = (
                f"{peptide_hit.getSequence().toString()}/"
                f"{peptide_hit.getCharge()}/{rt}/{spectrum_ref}"
            )

            return PSM(
                peptidoform=peptidoform,
                spectrum_id=spectrum_ref,
                run=IdXMLReader._get_run(protein_ids, peptide_id),
                is_decoy=is_decoy,
                score=peptide_hit.getScore(),
                precursor_mz=peptide_id.getMZ(),
                retention_time=rt,
                rank=peptide_hit.getRank() + 1,
                source="idXML",
                provenance_data={provenance_key: provenance_value},
            )
        except Exception as e:
            logging.error(f"Failed to parse PSM: {e}")
            return None

    def _parse_idxml(self) -> Tuple[List[oms.ProteinIdentification], List[oms.PeptideIdentification]]:
        """Parse the idXML file to extract protein and peptide identifications."""
        idxml_file = oms.IdXMLFile()
        proteins, peptides = [], []
        idxml_file.load(str(self.filename), proteins, peptides)
        return proteins, peptides

    def read_file(self) -> PSMList:
        """Read and parse all PSMs from the idXML file."""
        psm_list = []

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
        """Build a spectrum lookup index from an mzML file."""
        self._mzml_path = str(mzml_file) if isinstance(mzml_file, Path) else mzml_file
        self._exp, self._spec_lookup = OpenMSHelper.get_spectrum_lookup_indexer(self._mzml_path)
        logging.info(f"Built SpectrumLookup from {self._mzml_path}")

    def validate_psm_spectrum_references(self) -> SpectrumStats:
        """
        Validate PSM spectrum references and collect spectrum statistics.

        Returns:
            SpectrumStats object containing validation results

        Raises:
            ValueError: If spectrum lookup or PSMs are not initialized
        """
        if self._spec_lookup is None or self._exp is None or not self._psms:
            raise ValueError("Spectrum lookup or PSMs not initialized")

        stats = SpectrumStats()

        for psm in self._psms:
            spectrum = OpenMSHelper.get_spectrum_for_psm(psm, self._exp, self._spec_lookup)

            if spectrum is None:
                logging.error(f"Spectrum not found for PSM {psm}")
                stats.missing_spectra += 1
                continue

            peaks = spectrum.get_peaks()[0]
            if len(peaks) == 0:
                logging.warning(f"Empty spectrum found for PSM {psm}")
                stats.empty_spectra += 1

            ms_level = spectrum.getMSLevel()
            stats.ms_level_counts[ms_level] += 1

            if ms_level == 3:
                logging.info(f"MS level 3 spectrum found for PSM {psm}")

        if stats.missing_spectra or stats.empty_spectra:
            logging.error(
                f"Found {stats.missing_spectra} PSMs with missing spectra and "
                f"{stats.empty_spectra} PSMs with empty spectra"
            )

        logging.info(f"MS level distribution: {dict(stats.ms_level_counts)}")
        return stats