# Written by Jonas Scheid under the MIT license
# Contributions by Yasset Perez-Riverol and Dai Chengxin
# This script is part of the quantmsutils package

import logging
import os
from typing import List

import click
import pyopenms as oms
from ms2rescore import rescore
from psm_utils import PSMList
from psm_utils.io.idxml import IdXMLWriter

from quantmsrescore.annotator import Annotator

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

# Define scores that needs to be present in any peptide hit, if one PSM doesn't have these scores, it will be removed
# and Error will be printed. The scores are based on the search engine used to generate the idXML file
CORE_FEATURES = [
    # MS-GF+ Scores and E-values
    {"MS:1002050", "MS-GF:DeNovoScore"},  # MS-GF:DeNovoScore
    {"MS:1002049", "MS-GF:RawScore"},  # MS-GF:RawScore
    {"MS:1002052", "MS-GF:SpecEValue"},  # MS-GF:SpecEValue
    # Comet Scores and E-values
    {"MS:1002252", "Comet:xcorr"},  # Comet:Score
    {"MS:1002255", "Comet:spscore"},
    {"MS:1002253", "Comet:deltacn", "COMET:deltCn"},
    # Sage Scores and E-values
    {"ln(hyperscore)"},
    {"SAGE:matched_peaks"},
    {"SAGE:ln(delta_best)"},
    {"SAGE:ln(delta_next)"},
]

OPTIONAL_FEATURES = [
    # MSGF+ Scores and E-values
    {"MSGF:ScoreRatio"},
    {"MSGF:Energy"},
    {"MSGF:lnEValue"},
    {"MSGF:lnExplainedIonCurrentRatio"},
    {"MSGF:lnNTermIonCurrentRatio"},
    {"MSGF:lnCTermIonCurrentRatio"},
    {"MSGF:lnMS2IonCurrent"},
    {"MSGF:MeanErrorTop7"},
    {"MSGF:sqMeanErrorTop7"},
    {"MSGF:StdevErrorTop7"},
]


# def parse_cli_arguments_to_config(
#     config_file: str = None,
#     feature_generators: str = None,
#     ms2pip_model_dir: str = None,
#     ms2pip_model: str = None,
#     ms2_tolerance: float = None,
#     calibration_set_size: float = None,
#     rng: int = None,
#     processes: int = None,
#     spectrum_path: str = None,
#     lower_score_is_better: bool = None,
#     output_path: str = None,
#     log_level: str = None,
#     spectrum_id_pattern: str = None,
#     psm_id_pattern: str = None,
# ) -> dict:
#     if config_file is None:
#         config = json.load(importlib.resources.open_text(package_data, "config_default.json"))
#     else:
#         with open(config_file) as f:
#             config = json.load(f)
#     if feature_generators is not None:
#         feature_generators_list = feature_generators.split(",")
#         config["ms2rescore"]["feature_generators"] = {}
#         if "basic" in feature_generators_list:
#             config["ms2rescore"]["feature_generators"]["basic"] = {}
#         if "ms2pip" in feature_generators_list:
#             config["ms2rescore"]["feature_generators"]["ms2pip"] = {
#                 "model_dir": ms2pip_model_dir,
#                 "model": ms2pip_model,
#                 "ms2_tolerance": ms2_tolerance,
#             }
#         if "deeplc" in feature_generators_list:
#             config["ms2rescore"]["feature_generators"]["deeplc"] = {
#                 "deeplc_retrain": False,
#                 "calibration_set_size": calibration_set_size,
#             }
#         if "maxquant" in feature_generators_list:
#             config["ms2rescore"]["feature_generators"]["maxquant"] = {}
#         if "ionmob" in feature_generators:
#             config["ms2rescore"]["feature_generators"]["ionmob"] = {}
#
#     logging.info(
#         "Percolator or makapot rescoring happends outside of this tool, please check quantms workflow"
#     )
#
#     if ms2pip_model_dir is not None:
#         config["ms2rescore"]["ms2pip_model_dir"] = ms2pip_model_dir
#     if ms2pip_model is not None:
#         config["ms2rescore"]["ms2pip_model"] = ms2pip_model
#     if ms2_tolerance is not None:
#         config["ms2rescore"]["ms2_tolerance"] = ms2_tolerance
#     if calibration_set_size is not None:
#         config["ms2rescore"]["calibration_set_size"] = calibration_set_size
#     if rng is not None:
#         config["ms2rescore"]["rng"] = rng
#     if spectrum_path is not None:
#         config["ms2rescore"]["spectrum_path"] = spectrum_path
#     if id_decoy_pattern is not None:
#         config["ms2rescore"]["id_decoy_pattern"] = id_decoy_pattern
#     if lower_score_is_better is not None:
#         config["ms2rescore"]["lower_score_is_better"] = lower_score_is_better
#     if processes is None:
#         processes = 1  # Default to single process
#     config["ms2rescore"]["processes"] = processes
#     if output_path is not None:
#         config["ms2rescore"]["output_path"] = output_path
#     else:
#         raise ValueError("Output path must be specified.")
#     if log_level is not None:
#         config["ms2rescore"]["log_level"] = log_level
#     if spectrum_id_pattern is not None:
#         config["ms2rescore"]["spectrum_id_pattern"] = spectrum_id_pattern
#     if psm_id_pattern is not None:
#         config["ms2rescore"]["psm_id_pattern"] = psm_id_pattern
#
#     return config


def rescore_idxml(input_file, output_file, config) -> None:
    """Rescore PSMs in an idXML file and keep other information unchanged."""
    # Read PSMs
    reader = IDXMLReaderPatch(input_file)
    psm_list = reader.read_file()

    if reader.skip_invalid_psm != 0:
        logging.warning(f"Removed {reader.skip_invalid_psm} PSMs without search engine features!")
        # Synchronised acquisition of new peptide IDs after removing invalid PSMs
        peptide_ids = reader.new_peptide_ids
    else:
        peptide_ids = reader.oms_peptides

    # check if any spectrum is empty
    exp = oms.MSExperiment()
    oms.MzMLFile().load(config["ms2rescore"]["spectrum_path"], exp)
    empty_spectra = 0
    spec = []
    for spectrum in exp:
        peaks_tuple = spectrum.get_peaks()
        if len(peaks_tuple[0]) == 0 and spectrum.getMSLevel() == 2:
            logging.warning(f"{spectrum.getNativeID()} spectra don't have spectra information!")
            empty_spectra += 1
            continue
        spec.append(spectrum)

    if empty_spectra != 0:
        logging.warning(f"Removed {empty_spectra} spectra without spectra information!")
        exp.setSpectra(spec)
        output_dir = os.path.dirname(config["ms2rescore"]["output_path"])
        mzml_output = os.path.join(
            output_dir,
            os.path.splitext(os.path.basename(config["ms2rescore"]["spectrum_path"]))[0]
            + "_clear.mzML",
        )
        oms.MzMLFile().store(mzml_output, exp)
        config["ms2rescore"]["spectrum_path"] = mzml_output
        # TODO: Add cleanup of temporary file after processing

    # Rescore
    rescore(config, psm_list)

    # Filter out PeptideHits within PeptideIdentification(s) that could not be processed by all feature generators
    peptide_ids_filtered = filter_out_artifact_psms(psm_list, peptide_ids)

    # Write
    writer = IdXMLWriter(output_file, reader.oms_proteins, peptide_ids_filtered)
    writer.write_file(psm_list)


def filter_out_artifact_psms(
    psm_list: PSMList, peptide_ids: List[oms.PeptideIdentification]
) -> List[oms.PeptideIdentification]:
    """Filter out PeptideHits that could not be processed by all feature generators"""
    num_mandatory_features = max([len(psm.rescoring_features) for psm in psm_list])
    new_psm_list = PSMList(
        psm_list=[psm for psm in psm_list if len(psm.rescoring_features) == num_mandatory_features]
    )

    # get differing peptidoforms of both psm lists
    psm_list_peptides = set([next(iter(psm.provenance_data.items()))[1] for psm in psm_list])
    new_psm_list_peptides = set(
        [next(iter(psm.provenance_data.items()))[1] for psm in new_psm_list]
    )
    not_supported_peptides = psm_list_peptides - new_psm_list_peptides

    # no need to filter if all peptides are supported
    if len(not_supported_peptides) == 0:
        return peptide_ids
    # Create new peptide ids and filter out not supported peptides
    new_peptide_ids = []
    for peptide_id in peptide_ids:
        new_hits = []
        for hit in peptide_id.getHits():
            if hit.getSequence().toString() in not_supported_peptides:
                continue
            new_hits.append(hit)
        if len(new_hits) == 0:
            continue
        peptide_id.setHits(new_hits)
        new_peptide_ids.append(peptide_id)
    logging.info(
        f"Removed {len(psm_list_peptides) - len(new_psm_list_peptides)} PSMs. Peptides not supported: {not_supported_peptides}"
    )
    return new_peptide_ids


@click.command(
    "annotate",
    short_help="Annotate PSMs in an idXML file and keep other information unchanged.",
)
@click.option(
    "-p",
    "--psm_file",
    help="Path to PSM file (idXML)",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-s",
    "--spectrum_path",
    help="Path to MGF/mzML spectrum file or directory with spectrum files (default: derived from identification file)",
    required=True,
    type=click.Path(exists=True),
)
@click.option(
    "-o",
    "--output_path",
    help="Path and stem for output file names (default: derive from identification file)",
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
    default="",
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
@click.pass_context
def annotate(
    ctx,
    psm_file: str,
    spectrum_path,
    output_path: str,
    log_level,
    processes,
    feature_generators,
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

    Args
    ----
    ctx: Click context object.
    psm_file (str): Path to the input PSM file in idXML format.
    spectrum_path: Path to the spectrum file(s) in MGF/mzML format.
    output_path (str): Path and stem for output file names.
    log_level: Logging level for the operation.
    processes: Number of parallel processes to use.
    feature_generators: Comma-separated list of feature generators to use.
    ms2pip_model_dir: Directory path for the MS²PIP model.
    ms2pip_model: Name of the MS²PIP model to use.
    ms2_tolerance: Fragment mass tolerance in Daltons.
    calibration_set_size: Percentage size of the calibration set for DeepLC.
    id_decoy_pattern: Regular expression pattern for identifying decoys.
    lower_score_is_better: Boolean indicating if lower scores are better.
    config_file (str): Path to a configuration file.
    spectrum_id_pattern (str): Pattern for spectrum identification.
    psm_id_pattern (str): Pattern for PSM identification.
    """

    logging.getLogger().setLevel(log_level.upper())

    if output_path is None:
        output_path = psm_file.replace(".idXML", "_ms2rescore.idXML")

    annotator = Annotator(
        feature_generators=feature_generators,
        ms2pip_model=ms2pip_model,
        ms2pip_model_path=ms2pip_model_dir,
        ms2_tolerance=ms2_tolerance,
        calibration_set_size=calibration_set_size,
        processes=processes,
        id_decoy_pattern=id_decoy_pattern,
        lower_score_is_better=lower_score_is_better,
        log_level=log_level,
        spectrum_id_pattern=spectrum_id_pattern,
        psm_id_pattern=psm_id_pattern,
    )
    logging.info("MS²Rescore config:")
