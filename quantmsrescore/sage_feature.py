import click
import pandas as pd

from quantmsrescore.idxmlreader import IdXMLReader
from quantmsrescore.logging_config import get_logger, configure_logging
from quantmsrescore.openms import OpenMSHelper

# Configure logging with default settings
configure_logging()

# Get logger for this module
logger = get_logger(__name__)


@click.command("sage2feature")
@click.option("--idxml", "-i", help="Input idXML file")
@click.option("--output_file", "-o", help="Output idXML file")
@click.option("--feat_file", "-f", help="Input feature table file")
@click.pass_context
def add_sage_feature(ctx, idxml: str, output_file: str, feat_file: str):
    """
    Add extra features in features idXML. Adding extra feature in Sage isn't known input for PSMFeatureExtractor

    :param ctx: click context
    :param idxml: Original idXML file
    :param output_file: Outpuf file with the extra feature
    :param feat_file: Feature file from Sage
    :return: None
    """

    logger.info("Starting adding extra feature to {}".format(idxml))
    logger.info("Reading feature file")

    extra_feat = []
    feat = pd.read_csv(feat_file, sep="\t")
    for _, row in feat.iterrows():
        if row["feature_generator"] == "psm_file":
            continue
        else:
            extra_feat.append(row["feature_name"])

    logger.info("Reading idXML file")

    try:
        idxml_reader = IdXMLReader(idxml_filename=idxml)
        protein_ids = idxml_reader.oms_proteins
        peptide_ids = idxml_reader.oms_peptides
    except Exception as e:
        raise click.ClickException(f"Failed to read idXML file: {str(e)}")
    logger.info("Adding extra feature to idXML file")
    search_parameters = protein_ids[0].getSearchParameters()
    try:
        features = search_parameters.getMetaValue("extra_features")
    except Exception:
        features = ""
    extra_features = (features + "," if features else "") + ",".join(extra_feat)
    search_parameters.setMetaValue("extra_features", extra_features)
    protein_ids[0].setSearchParameters(search_parameters)

    logger.info("Writing idXML file")
    OpenMSHelper.write_idxml_file(
        filename=output_file, protein_ids=protein_ids, peptide_ids=peptide_ids
    )
    logger.info("Done")
