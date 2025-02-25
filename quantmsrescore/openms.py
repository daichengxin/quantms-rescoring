import logging
import re
from pathlib import Path
from typing import List, Union, Optional, Tuple

import numpy as np
import pyopenms as oms
from psm_utils import PSM
from pyopenms import PeptideIdentification, ProteinIdentification, SpectrumLookup, PeptideHit

OPENMS_DECOY_FIELD = "target_decoy"
SPECTRUM_PATTERN = r"(spectrum|scan)=(\d+)"


class OpenMSHelper:
    """
    This class should contain methods to help with OpenMS operations on PSM lists
    Features, parameters.
    """

    @staticmethod
    def count_decoys_targets(
        peptide_list: Union[List[PeptideIdentification], List[PeptideHit]],
    ) -> (int, int):
        """
        Count the number of decoy and target PSMs in the given list.

        This method iterates over a list of PSM objects, counting how many
        are labeled as 'target' or 'decoy' based on their rescoring features
        and the `is_decoy` attribute. It ensures that the counts from the
        rescoring features match the counts from the `is_decoy` attribute.

        Parameters
        ----------
        peptide_list (List[PSM]): A list of PeptideIdentification objects to be analyzed.

        Returns
        -------
        tuple: A tuple containing the count of decoy PSMs and the count of
        target PSMs.

        Raises
        -------
        ValueError: If the counts from the rescoring features do not match
        the counts from the `is_decoy` attribute.
        """

        openms_count_target = 0
        openms_count_decoy = 0

        for pep in peptide_list:
            if isinstance(pep, PeptideHit):
                if OpenMSHelper.is_decoy_peptide_hit(pep):
                    openms_count_decoy += 1
                else:
                    openms_count_target += 1
            else:
                for psm in pep.getHits():
                    if psm.metaValueExists(OPENMS_DECOY_FIELD):
                        if psm.getMetaValue(OPENMS_DECOY_FIELD) == "decoy":
                            openms_count_decoy += 1
                        else:
                            openms_count_target += 1

        if openms_count_decoy + openms_count_target == 0:
            logging.warning("No PSMs found; decoy percentage cannot be computed.")
            return 0, 0
        percentage_decoy = (openms_count_decoy / (openms_count_decoy + openms_count_target)) * 100
        logging.info(
            "Decoy percentage: %s, targets %s and decoys %s",
            percentage_decoy,
            openms_count_target,
            openms_count_decoy,
        )
        return openms_count_decoy, openms_count_target

    @staticmethod
    def get_psm_count(peptide_list: Union[List[PeptideIdentification], List[PeptideHit]]) -> int:
        """
        Count the number of PSMs in the given list.

        This method iterates over a list of PSM objects, counting the total
        number of PSMs.

        Parameters
        ----------
        peptide_list (List[PSM]): A list of PeptideIdentification objects to be analyzed.

        Returns
        -------
        int: The total number of PSMs in the list.
        """

        openms_count = 0

        for pep in peptide_list:
            if isinstance(pep, PeptideHit):
                openms_count += 1
            else:
                openms_count += len(pep.getHits())
        logging.info("Total PSMs: %s", openms_count)
        return openms_count

    @staticmethod
    def is_decoy_peptide_hit(peptide_hit: PeptideHit) -> bool:
        """
        Check if a PeptideHit is a decoy.

        This method checks if a PeptideHit is a decoy based on the
        'target_decoy' field in the PeptideHit.

        Parameters
        ----------
        peptide_hit (PeptideIdentification): A PeptideIdentification object to be checked.

        Returns
        -------
        bool: True if the PeptideHit is a decoy, False otherwise.
        """

        if peptide_hit.metaValueExists(OPENMS_DECOY_FIELD):
            return peptide_hit.getMetaValue(OPENMS_DECOY_FIELD) == "decoy"
        return False

    @staticmethod
    def get_spectrum_lookup_indexer(
        mzml_file: Union[str, Path],
    ) -> tuple[oms.MSExperiment, SpectrumLookup]:
        """
        Create a SpectrumLookup indexer from an mzML file.

        This method loads an mzML file into an MSExperiment object and
        initializes a SpectrumLookup object to read spectra using a
        specified regular expression pattern for scan numbers.

        Parameters
        ----------
        mzml_file : str
        The path to the mzML file to be loaded.

        Returns
        -------
        tuple: A tuple containing the MSExperiment object with the loaded
        """

        if isinstance(mzml_file, Path):
            mzml_file = str(mzml_file)

        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_file, exp)

        lookup = SpectrumLookup()
        lookup.readSpectra(exp, "scan=(?<SCAN>\\d+)")
        return exp, lookup

    @staticmethod
    def get_spectrum_for_psm(
        psm: Union[PSM, PeptideIdentification], exp: oms.MSExperiment, lookup: SpectrumLookup
    ) -> Union[None, oms.MSSpectrum]:
        spectrum_reference = ""
        if isinstance(psm, PSM):
            spectrum_reference = psm.spectrum_id
        elif isinstance(psm, PeptideIdentification):
            spectrum_reference = psm.getMetaValue("spectrum_reference")

        matches = re.findall(r"(spectrum|scan)=(\d+)", spectrum_reference)
        if not matches:
            psm_info = psm.provenance_data if hasattr(psm, "provenance_data") else "N/A"
            logging.warning(
                f"Missing or invalid spectrum reference for PSM {psm_info}, "
                f"skipping spectrum retrieval."
            )
            return None
        scan_number = int(matches[0][1])

        try:
            index = lookup.findByScanNumber(scan_number)
            spectrum = exp.getSpectrum(index)
            return spectrum
        except Exception as e:
            psm_info = psm.provenance_data if hasattr(psm, "provenance_data") else "N/A"
            logging.error(
                "Error while retrieving spectrum for PSM %s spectrum_ref %s: %s",
                psm_info,
                spectrum_reference,
                e,
            )
        return None

    @staticmethod
    def write_idxml_file(
        filename: Union[str, Path],
        peptide_ids: List[PeptideIdentification],
        protein_ids: List[ProteinIdentification],
    ) -> None:
        """
        Write protein and peptide identifications to an idXML file.

        Parameters
        ----------
        filename : Union[str, Path]
            The path to the idXML file to be written.
        peptide_ids : List[PeptideIdentification]
            A list of PeptideIdentification objects to be written to the file.
        protein_ids : List[ProteinIdentification]
            A list of ProteinIdentification objects to be written to the file.

        """

        if isinstance(filename, Path):
            filename = str(filename)

        id_data = oms.IdXMLFile()
        id_data.store(filename, protein_ids, peptide_ids)

    @staticmethod
    def get_peaks_by_scan(
        scan_number: int, exp: oms.MSExperiment, lookup: SpectrumLookup
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """
        Get spectrum data for a given scan number

        Parameters
        ----------
            scan_number: The scan number to look up
            exp: The MSExperiment object containing the spectra
            lookup: The SpectrumLookup object to find the spectrum

        Returns
        -------
            Tuple of (mz_array, intensity_array) if found, None if not found
        """
        try:
            index = lookup.findByScanNumber(scan_number)
            spectrum = exp.getSpectrum(index)
            return spectrum.get_peaks()
        except IndexError:
            logging.warning(f"Scan number {scan_number} not found")
            return None

    @staticmethod
    def get_ms_level(
        psm_hit: PeptideIdentification, spec_lookup: oms.SpectrumLookup, exp: oms.MSExperiment
    ) -> int:
        spectrum = OpenMSHelper.get_spectrum_for_psm(psm_hit, exp, spec_lookup)
        if spectrum is None:
            return -1
        return spectrum.getMSLevel()

    @staticmethod
    def get_psm_hash_unique_id(peptide_hit: PeptideIdentification, psm_hit: PeptideHit) -> str:
        """
        Generate a unique hash identifier for a PSM.

        This method constructs a unique hash string for a given PSM by
        combining the peptide sequence, charge, retention time, and
        spectrum reference.

        Parameters
        ----------
        peptide_hit : PeptideIdentification
            The PeptideIdentification object containing metadata for the PSM.
        psm_hit : PeptideHit
            The PeptideHit object representing the PSM.
        spectrum_file : str
            The path to the spectrum file containing the PSM. This is needed in case multiple files are
            provided in the same run.

        Returns
        -------
        str
            A unique hash string for the PSM.
        """

        spectrum_ref = peptide_hit.getMetaValue("spectrum_reference")
        rank = psm_hit.getRank()
        rt = peptide_hit.getRT()
        sequence = psm_hit.getSequence().toString()
        charge = psm_hit.getCharge()
        unique_hash = f"{spectrum_ref}_{sequence}_{rt}_{charge}_{rank}"
        return unique_hash

    @staticmethod
    def get_str_metavalue_round(metavalue: float):
        """
        Get a string representation of a metadata value, rounded to 4 decimal places.

        Parameters
        ----------
        metavalue : float
            The metadata value to be converted to a string.

        Returns
        -------
        str
            A string representation of the metadata value, rounded to 4 decimal places.
        """
        if np.isnan(metavalue) or np.isinf(metavalue):
            return "0.0"
        return "{:.4f}".format(metavalue)
