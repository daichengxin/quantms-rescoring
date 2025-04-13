import re
from dataclasses import dataclass

import click
import numpy as np
from scipy.stats import entropy

from quantmsrescore.idxmlreader import IdXMLReader
from quantmsrescore.logging_config import get_logger, configure_logging
from quantmsrescore.openms import OpenMSHelper

# Configure logging with default settings
configure_logging()

# Get logger for this module
logger = get_logger(__name__)


@dataclass
class SpectrumMetrics:
    """Data class to hold spectrum analysis metrics"""

    snr: float
    spectral_entropy: float
    fraction_tic_top_10: float
    weighted_std_mz: float

    def as_dict(self) -> dict:
        """Convert metrics to a dictionary with proper prefixes"""

        return {
            "Quantms:Snr": OpenMSHelper.get_str_metavalue_round(self.snr),
            "Quantms:SpectralEntropy": OpenMSHelper.get_str_metavalue_round(self.spectral_entropy),
            "Quantms:FracTICinTop10Peaks": OpenMSHelper.get_str_metavalue_round(
                self.fraction_tic_top_10
            ),
            "Quantms:WeightedStdMz": OpenMSHelper.get_str_metavalue_round(self.weighted_std_mz),
        }


class SpectrumAnalyzer:
    """Class to handle spectrum analysis operations"""

    @staticmethod
    def compute_signal_to_noise(intensities: np.ndarray) -> float:
        """
        Compute the signal-to-noise ratio for a given array of intensities.

        Parameters
        ----------
        intensities : np.ndarray
            Array of intensity values.

        Returns
        -------
        float
            The signal-to-noise ratio calculated as the maximum intensity
            divided by the root mean square deviation of the intensities.

        Raises
        ------
        ValueError
            If the input intensity array is empty.
        """
        if len(intensities) == 0:
            raise ValueError("Empty intensity array provided")

        rmsd = np.sqrt(np.mean(np.square(intensities)))
        if rmsd == 0:
            return 0

        return np.max(intensities) / rmsd

    @staticmethod
    def compute_spectrum_metrics(
        mz_array: np.ndarray, intensity_array: np.ndarray
    ) -> SpectrumMetrics:
        """
        Compute various spectrum metrics including signal-to-noise ratio,
        spectral entropy, fraction of total ion current in the top 10 peaks,
        and weighted standard deviation of m/z values.

        Parameters
        ----------
            mz_array (np.ndarray): Array of m/z values.
            intensity_array (np.ndarray): Array of intensity values.

        Returns
        -------
            SpectrumMetrics: An instance containing computed spectrum metrics.

        Raises
        ------
            ValueError: If input arrays are empty, have different lengths, or
                        if the total ion current is zero.
        """
        if len(intensity_array) == 0 or len(mz_array) == 0:
            raise ValueError("Empty arrays provided")

        if len(mz_array) != len(intensity_array):
            raise ValueError("mz_array and intensity_array must have same length")

        # Total Ion Current
        tic = np.sum(intensity_array)
        if tic == 0:
            raise ValueError("Total ion current is zero")

        # Normalized intensities
        normalized_intensities = intensity_array / tic

        # Calculate all metrics
        snr = SpectrumAnalyzer.compute_signal_to_noise(intensity_array)
        spectral_entropy = entropy(normalized_intensities)

        # Top 10 peaks analysis
        top_n_peaks = np.sort(intensity_array)[-10:]
        fraction_tic_top_10 = np.sum(top_n_peaks) / tic

        # Weighted m/z calculations
        weighted_mean_mz = np.sum(mz_array * normalized_intensities)
        weighted_std_mz = np.sqrt(
            np.sum(normalized_intensities * (mz_array - weighted_mean_mz) ** 2)
        )

        return SpectrumMetrics(snr, spectral_entropy, fraction_tic_top_10, weighted_std_mz)


@click.command("spectrum2feature")
@click.option(
    "--mzml",
    type=click.Path(exists=True),
    required=True,
    help="Path to the mass spectrometry file",
)
@click.option(
    "--idxml",
    type=click.Path(exists=True),
    required=True,
    help="Path to the idXML file with PSMs corresponding to the mzML file",
)
@click.option(
    "--output",
    type=click.Path(),
    required=True,
    help="Path for the output idXML file with computed metrics",
)
@click.pass_context
def spectrum2feature(ctx, mzml: str, idxml: str, output: str) -> None:
    """
    Command-line tool to compute spectrum metrics and update idXML files.

    This command processes an mzML file and an idXML file containing PSMs,
    computes spectrum metrics for each peptide identification, and updates
    the idXML file with these metrics.

    Parameters
    ----------
    ctx : click.Context
        The Click context object.
    mzml : str
        Path to the mzML file containing mass spectrometry deeplc_models.
    idxml : str
        Path to the idXML file with PSMs corresponding to the mzML file.
    output : str
        Path for the output idXML file with computed metrics.

    Raises
    ------
    ValueError
        If no protein identifications are found in the idXML file.
    """
    logger.info(f"Processing mzML file: {mzml}")

    idxml_reader = IdXMLReader(idxml_filename=idxml)
    idxml_reader.build_spectrum_lookup(mzml, check_unix_compatibility=True)
    protein_ids = idxml_reader.oms_proteins
    peptide_ids = idxml_reader.oms_peptides

    if not protein_ids:
        raise ValueError("No protein identifications found in idXML file")

    result_peptides = []
    for peptide in peptide_ids:
        spectrum_reference = peptide.getMetaValue("spectrum_reference")
        scan_matches = re.findall(r"(spectrum|scan)=(\d+)", spectrum_reference)

        if not scan_matches:
            logger.warning(f"Could not parse scan number from reference: {spectrum_reference}")
            continue

        scan_number = int(scan_matches[0][1])
        spectrum_data = OpenMSHelper.get_peaks_by_scan(scan_number, idxml_reader.exp, idxml_reader.spec_lookup)

        if spectrum_data is None:
            continue

        try:
            mz_array, intensity_array = spectrum_data
            metrics = SpectrumAnalyzer.compute_spectrum_metrics(
                np.array(mz_array), np.array(intensity_array)
            )

            # Update peptide hits with metrics
            for hit in peptide.getHits():
                for key, value in metrics.as_dict().items():
                    hit.setMetaValue(key, value)
                peptide.setHits([hit])

            result_peptides.append(peptide)

        except ValueError as e:
            logger.error(f"Error processing scan {scan_number}: {str(e)}")
            continue

    # Update search parameters with new features
    search_parameters = protein_ids[0].getSearchParameters()
    existing_features = search_parameters.getMetaValue("extra_features")
    new_features = ",".join(metrics.as_dict().keys())
    extra_features = f"{existing_features},{new_features}" if existing_features else new_features
    search_parameters.setMetaValue("extra_features", extra_features)
    protein_ids[0].setSearchParameters(search_parameters)

    # Save results
    OpenMSHelper.write_idxml_file(
        filename=output, protein_ids=protein_ids, peptide_ids=result_peptides
    )
    logger.info(f"Results saved to: {output}")
