# Get logger for this module
from quantmsrescore.logging_config import get_logger
from collections import defaultdict
from pathlib import Path
from typing import Union, List, Optional, Dict, Tuple, DefaultDict
from quantmsrescore.exceptions import MS3NotSupportedException
import pyopenms as oms
from pyopenms import IDFilter

from quantmsrescore.openms import OpenMSHelper

logger = get_logger(__name__)


class SpectrumStats:
    """Statistics about spectrum analysis."""

    def __init__(self):
        self.missing_spectra: int = 0
        self.empty_spectra: int = 0
        self.invalid_score: int = 0
        self.ms_level_counts: DefaultDict[int, int] = defaultdict(int)
        self.ms_level_dissociation_method: Dict[Tuple[int, str], int] = {}


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

    def __init__(self, idxml_filename: Union[Path, str]) -> None:
        """
        Initialize IdXMLReader with the specified idXML file.

        Parameters
        ----------
        idxml_filename : Union[Path, str]
            Path to the idXML file to be read and parsed.
        """
        self.filename = Path(idxml_filename)
        self.oms_proteins, self.oms_peptides = self._parse_idxml()
        self.spec_lookup = None
        self.exp = None

        # Private properties for spectrum lookup
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
        self.exp, self.spec_lookup = OpenMSHelper.get_spectrum_lookup_indexer(self._mzml_path)
        logger.info(f"Built SpectrumLookup from {self._mzml_path}")

    def psm_clean(self, remove_missing_spectrum: bool = True, only_ms2: bool = True) -> SpectrumStats:
        """
        Validate spectrum references for peptide identifications and filter based on criteria.

        This method validates each peptide identification by checking if its referenced
        spectrum exists and has peaks. It also tracks MS level statistics, dissociation
        methods and removes peptide identifications with missing/empty spectra or
        those that are not MS2 level and invalid score.

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

        if self.spec_lookup is None or self.exp is None:
            raise ValueError("Spectrum lookup or PSMs not initialized")

        self._stats = SpectrumStats()

        new_peptide_ids = []
        peptide_removed = 0
        search_engine = self.oms_proteins[0].getSearchEngine()

        for peptide_id in self.oms_peptides:
            spectrum = OpenMSHelper.get_spectrum_for_psm(peptide_id, self.exp, self.spec_lookup)
            spectrum_reference = OpenMSHelper.get_spectrum_reference(peptide_id)

            missing_spectrum, empty_spectrum, invalid_score = False, False, False
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

            # Removed the Hit when Sage output Inf value for poisson score
            if search_engine == "Sage":
                hits_number = len(peptide_id.getHits())
                new_hits = []
                for hit in peptide_id.getHits():
                    score = hit.getMetaValue("SAGE:ln(-poisson)")
                    if score == "inf":
                        hits_number -= 1
                        logger.warning(f"Invalid PSM score found for PSM {spectrum_reference}")
                    else:
                        new_hits.append(hit)
                if hits_number == 0:
                    invalid_score = True
                    self._stats.invalid_score += 1
                elif hits_number < len(peptide_id.getHits()):
                    peptide_id.setHits(new_hits)

            if (remove_missing_spectrum and (missing_spectrum or empty_spectrum or invalid_score)) or (
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

        self._log_spectrum_statistics()

        if only_ms2 and self._stats.ms_level_counts.get(3, 0) > 0:
            ms2_dissociation_methods = self._stats.ms_level_dissociation_method.get((2, "HCD"), 0)
            logger.error(
                "MS3 spectra found in MS2-only mode, please filter your search for MS2 or dissociation method: {}".format(
                    ms2_dissociation_methods
                )
            )
            raise MS3NotSupportedException("MS3 spectra found in MS2-only mode")

        return self._stats

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

    def write_idxml_file(self, filename: Union[str, Path]) -> None:
        """
        Write processed data to idXML file.

        Parameters
        ----------
        filename : Union[str, Path]
            Path where the processed idXML file will be written.

        Raises
        ------
        Exception
            If writing the file fails.
        """
        try:
            out_path = Path(filename)
            OpenMSHelper.write_idxml_file(
                filename=out_path,
                protein_ids=self.openms_proteins,
                peptide_ids=self.openms_peptides,
            )
            logger.info(f"Processed idXML file written to {out_path}")
        except Exception as e:
            logger.error(f"Failed to write Processed idXML file: {str(e)}")
            raise
