import logging
from collections import defaultdict
from pathlib import Path
from typing import Union, Iterable, List, Optional, Dict, Tuple, DefaultDict, Set
from warnings import filterwarnings

import pyopenms as oms
from psm_utils import PSM, PSMList
from psm_utils.io.idxml import IdXMLReader

from quantmsrescore.constants import OPENMS_DISSOCIATION_METHODS_PATCH
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
    ms_level_dissociation_method = {}


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

    def __init__(
        self, idexml_filename: Union[Path, str], mzml_file: Union[str, Path] = None
    ) -> None:
        """
        Initialize the IdXMLRescoringReader with the specified idXML file.

        Parameters
        ----------
        idexml_filename : Union[Path, str]
            The path to the idXML file to be processed.

        Attributes
        ----------
        filename : Path
            The path to the idXML file.
        high_score_better : Optional[bool]
            Indicates if higher scores are better.
        skip_invalid_psm : int
            Counter for skipped invalid PSMs.
        new_peptide_ids : List[oms.PeptideIdentification]
            List to store new peptide identifications.
        _mzml_path : Optional[str]
            Path to the mzML file, if specified.
        _psms : Optional[PSMList]
            List of parsed PSMs.
        _spec_lookup : Optional[oms.SpectrumLookup]
            Spectrum lookup index.
        _exp : Optional[oms.MSExperiment]
            MS experiment data.
        oms_proteins : List[oms.ProteinIdentification]
            List of protein identifications from the idXML file.
        oms_peptides : List[oms.PeptideIdentification]
            List of peptide identifications from the idXML file.
        """
        self.filename = Path(idexml_filename)
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

        # If mzML file is provided, build SpectrumLookup
        if mzml_file:
            self.build_spectrum_lookup(mzml_file)

    @property
    def spectrum_path(self) -> Union[str, Path]:
        return self._mzml_path

    @property
    def psms(self) -> Optional[PSMList]:
        return self._psms

    @psms.setter
    def psms(self, psm_list: PSMList) -> None:
        if not isinstance(psm_list, PSMList):
            raise TypeError("psm_list must be an instance of PSMList")
        self._psms = psm_list

    @property
    def openms_proteins(self) -> List[oms.ProteinIdentification]:
        return self.oms_proteins

    @property
    def openms_peptides(self) -> List[oms.PeptideIdentification]:
        return self.oms_peptides

    def __iter__(self) -> Iterable[PSM]:
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
        Analyze the coverage of scores across peptide hits.

        This method calculates the total number of hits for each score
        present in the peptide identifications and determines the number
        of missing scores by comparing against the total number of hits.

        Returns
        -------
        Dict[str, ScoreStats]
            A dictionary mapping score names to their respective statistics,
            including total hits and missing counts.
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

    def _log_score_coverage(self, score_stats: Dict[str, ScoreStats]) -> None:
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
        """
        Parse a peptide-spectrum match (PSM) from given protein and peptide data.

        This static method extracts relevant information from the provided
        protein, peptide, and peptide hit objects to construct a PSM object.
        It handles the creation of provenance tracking data and logs errors
        if parsing fails.

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
            peptidoform = IdXMLReader._parse_peptidoform(
                peptide_hit.getSequence().toString(), peptide_hit.getCharge()
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

    def _parse_idxml(
        self,
    ) -> Tuple[List[oms.ProteinIdentification], List[oms.PeptideIdentification]]:
        """
        Parse the idXML file to extract protein and peptide identifications.

        Returns
        -------
            Tuple[List[oms.ProteinIdentification], List[oms.PeptideIdentification]]:
            A tuple containing lists of protein and peptide identifications.
        """
        idxml_file = oms.IdXMLFile()
        proteins, peptides = [], []
        idxml_file.load(str(self.filename), proteins, peptides)
        return proteins, peptides

    def build_psm_index(self, only_ms2: bool = True) -> PSMList:
        """
        Read and parse the idXML file to extract PSMs.

        This method processes peptide identifications from the idXML file,
        determines the score direction if not already set, and parses each
        peptide hit into a PSM object. It logs warnings for inconsistent
        score directions and returns a list of valid PSMs.

        Parameters
        ----------
        only_ms2 : bool, optional
            Flag to filter for MS2 spectra only, by default True.

        Returns
        -------
            PSMList: A list of parsed PSM objects.
        """
        psm_list = []

        if only_ms2 and self._spec_lookup is None:
            logging.warning("Spectrum lookup not initialized, cannot filter for MS2 spectra")
            only_ms2 = False

        only_ms2 = only_ms2 and self._spec_lookup is not None

        for peptide_id in self.oms_peptides:
            if self.high_score_better is None:
                self.high_score_better = peptide_id.isHigherScoreBetter()
            elif self.high_score_better != peptide_id.isHigherScoreBetter():
                logging.warning("Inconsistent score direction found in idXML file")

            for psm_hit in peptide_id.getHits():
                if (
                    only_ms2
                    and OpenMSHelper.get_ms_level(psm_hit, self._spec_lookup, self._exp) != 2
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
        logging.info(f"Loaded {len(self._psms)} PSMs from {self.filename}")
        return self._psms

    def build_spectrum_lookup(self, mzml_file: Union[str, Path]) -> None:
        """
        Build a SpectrumLookup indexer from an mzML file.

        This method initializes the spectrum lookup and MS experiment data
        from the specified mzML file, storing them as class attributes for
        further processing and analysis.

        Parameters
        ----------
        mzml_file : Union[str, Path]
            The path to the mzML file to be processed.
        """
        self._mzml_path = str(mzml_file) if isinstance(mzml_file, Path) else mzml_file
        self._exp, self._spec_lookup = OpenMSHelper.get_spectrum_lookup_indexer(self._mzml_path)
        logging.info(f"Built SpectrumLookup from {self._mzml_path}")

    def validate_psm_spectrum_references(
        self, remove_missing_spectrum: bool = True, only_ms2: bool = True
    ) -> SpectrumStats:
        """
        Validate and update PSM spectrum references.

        This method checks the spectrum references for each PSM in the list,
        updating statistics on missing and empty spectra, MS level counts, and
        dissociation methods. It optionally removes PSMs with missing or empty
        spectra and filters for MS2 spectra only.

        Parameters
        ----------
        remove_missing_spectrum : bool, optional
            If True, remove PSMs with missing or empty spectra, by default True.
        only_ms2 : bool, optional
            If True, filter for MS2 spectra only, by default True.

        Returns
        -------
        SpectrumStats
            An object containing statistics about the spectrum validation process.

        Raises
        ------
        ValueError
            If the spectrum lookup or PSMs are not initialized.
        """
        if self._spec_lookup is None or self._exp is None or not self._psms:
            raise ValueError("Spectrum lookup or PSMs not initialized")

        stats = SpectrumStats()

        for psm in self._psms[:]:  # Iterate over a copy to allow safe removal
            spectrum = OpenMSHelper.get_spectrum_for_psm(psm, self._exp, self._spec_lookup)
            missing_spectra, empty_spectra = False, False
            ms_level = 2

            if spectrum is None:
                logging.error(f"Spectrum not found for PSM {psm}")
                stats.missing_spectra += 1
                missing_spectra = True
            else:
                peaks = spectrum.get_peaks()[0]
                if not peaks:
                    logging.warning(f"Empty spectrum found for PSM {psm}")
                    empty_spectra = True
                    stats.empty_spectra += 1

                ms_level = spectrum.getMSLevel()
                stats.ms_level_counts[ms_level] += 1

                for precursor in spectrum.getPrecursors():
                    for method_index in precursor.getActivationMethods():
                        if 0 <= method_index < len(OPENMS_DISSOCIATION_METHODS_PATCH):
                            method = (
                                ms_level,
                                list(OPENMS_DISSOCIATION_METHODS_PATCH[method_index].keys())[0],
                            )
                            stats.ms_level_dissociation_method[method] = (
                                stats.ms_level_dissociation_method.get(method, 0) + 1
                            )
                        else:
                            logging.warning(f"Unknown dissociation method index {method_index}")

                if ms_level == 3:
                    logging.info(
                        f"MS level 3 spectrum found for PSM {psm}, please be aware, "
                        "ms2pip models are not trained on MS3 spectra"
                    )

            if (remove_missing_spectrum and (missing_spectra or empty_spectra)) or (
                only_ms2 and ms_level != 2
            ):
                self._psms.remove(psm)

        if stats.missing_spectra or stats.empty_spectra:
            logging.error(
                f"Found {stats.missing_spectra} PSMs with missing spectra and "
                f"{stats.empty_spectra} PSMs with empty spectra"
            )

        logging.info(f"MS level distribution: {dict(stats.ms_level_counts)}")
        print("Dissociation Method Distribution", stats.ms_level_dissociation_method)
        return stats
