# Get logger for this module
from quantmsrescore.logging_config import get_logger

logger = get_logger(__name__)

from warnings import filterwarnings

filterwarnings(
    "ignore",
    message="OPENMS_DATA_PATH environment variable already exists",
    category=UserWarning,
    module="pyopenms",
)



import click

from quantmsrescore.logging_config import configure_logging
from quantmsrescore.utils import IdXMLReader

# Configure logging with default settings
configure_logging()


@click.command(
    "psm_feature_clean",
    short_help="Annotate PSMs in an idXML file using ms2rescore features.",
)
@click.option(
    "-i",
    "--idxml",
    help="Path to the idxml containing the PSMs from OpenMS",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-s",
    "--mzml",
    help="Path to the mzML file containing the spectra use for identification",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--output",
    help="Path the output idxml for the processed PSMs",
)
@click.pass_context
def psm_feature_clean(
    ctx,
    idxml: str,
    mzml: str,
    output: str,
):
    """
    Annotate PSMs in an idXML file with additional features using specified models.

    This command-line interface (CLI) command processes a PSM file by remove invalid PSMs

    Parameters
    ----------

    ctx : click.Context
        The Click context object.
    idxml : str
        Path to the idXML file containing the PSMs.
    mzml : str
        Path to the mzML file containing the mass spectrometry.
    output : str
        Path to the output idXML file with processed PSMs.
    """

    id_reader = IdXMLReader(idxml)
    id_reader.build_spectrum_lookup(mzml, check_unix_compatibility=True)
    id_reader.psm_clean()
    id_reader.write_idxml_file(output)

