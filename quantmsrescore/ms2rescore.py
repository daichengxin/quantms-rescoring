# Written by Jonas Scheid under the MIT license
# Contributions by Yasset Perez-Riverol and Dai Chengxin
# This script is part of the quantmsutils package

import logging
import click

from quantmsrescore.annotator import Annotator

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")


@click.command(
    "msrescore2feature",
    short_help="Annotate PSMs in an idXML file using ms2rescore features.",
)
@click.option(
    "-p",
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
    "--output_path",
    help="Path the output idxml for the annotated PSMs",
)
@click.option("-l", "--log_level", help="Logging level (default: `info`)", default="info")
@click.option(
    "-n",
    "--processes",
    help="Number of parallel processes available to MS²Rescore",
    type=int,
    default=16,
)
@click.option(
    "-fg",
    "--feature_generators",
    help="Comma-separated list of feature generators to use (default: `ms2pip,deeplc`). See rescoring doc for further information",
    default="ms2pip,deeplc",
)
@click.option(
    "--only_features",
    help="Comma-separated list of features to use for annotation (read docs for default)",
)
@click.option(
    "-pipm",
    "--ms2pip_model",
    help="MS²PIP model (default: `HCD2021`)",
    type=str,
    default="HCD2021",
)
@click.option(
    "-md",
    "--ms2pip_model_dir",
    help="The path of MS²PIP model (default: `./`)",
    type=str,
    default="./",
)
@click.option(
    "-ms2tol",
    "--ms2_tolerance",
    help="Fragment mass tolerance [Da](default: `0.05`)",
    type=float,
    default=0.05,
)
@click.option(
    "-cs",
    "--calibration_set_size",
    help="Percentage of number of calibration set for DeepLC (default: `0.15`)",
    default=0.15,
)
@click.option(
    "-d",
    "--id_decoy_pattern",
    help="Regex decoy pattern (default: `DECOY_`)",
    default="^DECOY_",
)
@click.option(
    "--lower_score_is_better",
    is_flag=True,
    help="Specify if lower scores are better",
    default=True,
)
@click.option(
    "--spectrum_id_pattern",
    help="Pattern for spectrum identification",
    type=str,
    default="(.*)",
)
@click.option(
    "--psm_id_pattern",
    help="Pattern for PSM identification",
    type=str,
    default="(.*)",
)
@click.pass_context
def msrescore2feature(
    ctx,
    idxml: str,
    mzml,
    output_path: str,
    log_level,
    processes,
    feature_generators,
    only_features,
    ms2pip_model_dir,
    ms2pip_model,
    ms2_tolerance,
    calibration_set_size,
    id_decoy_pattern,
    lower_score_is_better,
    spectrum_id_pattern: str,
    psm_id_pattern: str,
):
    """
    Annotate PSMs in an idXML file with additional features using specified models.

    This command-line interface (CLI) command processes a PSM file by adding
    annotations from the MS²PIP and DeepLC models, among others, while preserving
    existing information. It supports various options for specifying input and
    output paths, logging levels, and feature generation configurations.

    Parameters
    ----------

    ctx : click.Context
        The Click context object.
    idxml : str
        Path to the idXML file containing the PSMs.
    mzml : str
        Path to the mzML file containing the mass spectrometry deeplc_models.
    output_path : str
        Path to the output idXML file with annotated PSMs.
    log_level : str
        The logging level for the CLI command.
    processes : int
        The number of parallel processes available for MS²Rescore.
    feature_generators : str
        Comma-separated list of feature generators to use for annotation.
    only_features : str
        Comma-separated list of features to use for annotation.
    ms2pip_model_dir : str
        Path to the directory containing the MS²PIP models.
    ms2pip_model : str
        The MS²PIP model to use for annotation.
    ms2_tolerance : float
        The tolerance for MS²PIP annotation.
    calibration_set_size : float
        The fraction of deeplc_models used for calibration.
    id_decoy_pattern : str
        The regex pattern to identify decoy IDs.
    lower_score_is_better : bool
        Whether a lower score indicates a better match.
    spectrum_id_pattern : str
        The regex pattern for spectrum IDs.
    psm_id_pattern : str
        The regex pattern for PSM IDs.
    """

    annotator = Annotator(
        feature_generators=feature_generators,
        only_features=only_features,
        ms2pip_model=ms2pip_model,
        ms2pip_model_path=ms2pip_model_dir,
        ms2_tolerance=ms2_tolerance,
        calibration_set_size=calibration_set_size,
        skip_deeplc_retrain=True,
        processes=processes,
        id_decoy_pattern=id_decoy_pattern,
        lower_score_is_better=lower_score_is_better,
        log_level=log_level.upper(),
        spectrum_id_pattern=spectrum_id_pattern,
        psm_id_pattern=psm_id_pattern,
    )
    annotator.build_idxml_data(idxml, mzml)
    annotator.annotate()

    if output_path:
        annotator.write_idxml_file(output_path)
