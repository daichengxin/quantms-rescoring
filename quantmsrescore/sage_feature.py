import logging
from warnings import filterwarnings

import click
import pandas as pd

from quantmsrescore.idxmlreader import IdXMLReader
from quantmsrescore.openms import OpenMSHelper

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# Suppress OpenMS warning about a data path
filterwarnings(
    "ignore",
    message="OPENMS_DATA_PATH environment variable already exists",
    category=UserWarning,
    module="pyopenms",
)


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

    logging.info("Starting adding extra feature to {}".format(idxml))
    logging.info("Reading feature file")

    extra_feat = []
    feat = pd.read_csv(feat_file, sep="\t")
    for _, row in feat.iterrows():
        if row["feature_generator"] == "psm_file":
            continue
        else:
            extra_feat.append(row["feature_name"])

    logging.info("Reading idXML file")

    idxml_reader = IdXMLReader(idexml_filename=idxml)
    protein_ids = idxml_reader.oms_proteins
    peptide_ids = idxml_reader.oms_peptides

    logging.info("Adding extra feature to idXML file")
    search_parameters = protein_ids[0].getSearchParameters()
    features = search_parameters.getMetaValue("extra_features")
    extra_features = features + "," + ",".join(extra_feat)
    search_parameters.setMetaValue("extra_features", extra_features)
    protein_ids[0].setSearchParameters(search_parameters)

    logging.info("Writing idXML file")
    OpenMSHelper.write_idxml_file(
        filename=output_file, protein_ids=protein_ids, peptide_ids=peptide_ids
    )
    logging.info("Done")
