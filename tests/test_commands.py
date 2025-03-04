import logging
from pathlib import Path

import pytest
from click.testing import CliRunner

from quantmsrescore.annotator import Annotator
from quantmsrescore.idxmlreader import IdXMLRescoringReader
from quantmsrescore.openms import OpenMSHelper
from quantmsrescore.rescoring import cli

TESTS_DIR = Path(__file__).parent

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@pytest.mark.skip(reason="This is for local test in big datasets, kipping for now")
def test_ms2rescore():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "ms2rescore",
            "--psm_file",
            "{}/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML".format(
                TESTS_DIR
            ),
            "--spectrum_path",
            "{}/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML".format(
                TESTS_DIR
            ),
            "--processes",
            "2",
            "--ms2pip_model",
            "HCD2021",
            "--feature_generators",
            "'ms2pip,deeplc'",
            "--id_decoy_pattern",
            "^rev",
        ],
    )
    assert result.exit_code == 0


@pytest.mark.skip(reason="This is for local test in big datasets, kipping for now")
def test_ms2rescore_failing():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "ms2rescore",
            "--psm_file",
            "{}/test_data/dae1cb16fb57893b94bfcb731b2bf7/Instrument1_sample14_S1R10_042116_Fr12_msgf.idXML".format(
                TESTS_DIR
            ),
            "--spectrum_path",
            "{}/test_data/dae1cb16fb57893b94bfcb731b2bf7/Instrument1_sample14_S1R10_042116_Fr12.mzML".format(
                TESTS_DIR
            ),
            "--processes",
            "2",
            "--ms2pip_model",
            "TMT",
            "--id_decoy_pattern",
            "^DECOY_",
            "--ms2_tolerance",
            "0.4",
            "--calibration_set_size",
            "0.15",
            "--feature_generators",
            "deeplc,ms2pip",
        ],
    )
    assert result.exit_code == 0


def test_idxmlreader():

    idxml_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML"
    )

    mzml_file = (
        TESTS_DIR / "test_data" / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML"
    )

    idxml_reader = IdXMLRescoringReader(idexml_filename=idxml_file, mzml_file=mzml_file)
    logging.info("Loaded %s PSMs from %s", len(idxml_reader.psms), idxml_file)
    assert len(idxml_reader.psms) == 5346

    openms_helper = OpenMSHelper()
    decoys, targets = openms_helper.count_decoys_targets(idxml_reader.oms_peptides)
    logging.info(
        "Loaded %s PSMs from %s, %s decoys and %s targets",
        len(idxml_reader.psms),
        idxml_file,
        decoys,
        targets,
    )
    assert decoys == 1923
    assert targets == 3423
    stats = idxml_reader.stats
    assert stats.missing_spectra == 4


def test_annotator_train_rt():

    idxml_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML"
    )

    mzml_file = (
        TESTS_DIR / "test_data" / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML"
    )

    annotator = Annotator(
        feature_generators="ms2pip,deeplc",
        ms2pip_model="TMT",
        ms2pip_model_path="models",
        ms2_tolerance=0.05,
        calibration_set_size=0.15,
        skip_deeplc_retrain=False,
        processes=2,
        id_decoy_pattern="^DECOY_",
        lower_score_is_better=True,
        log_level="INFO",
        spectrum_id_pattern="(.*)",
        psm_id_pattern="(.*)",
    )
    annotator.build_idxml_data(idxml_file, mzml_file)
    annotator.annotate()

    output_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet_rescored.idXML"
    )

    annotator.write_idxml_file(output_file)


def test_idxmlreader_filtering():

    idxml_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML"
    )

    mzml_file = (
        TESTS_DIR / "test_data" / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML"
    )

    annotator = Annotator(
        feature_generators="ms2pip,deeplc",
        only_features="DeepLC:RtDiff,DeepLC:PredictedRetentionTimeBest,Ms2pip:DotProd",
        ms2pip_model="TMT",
        ms2pip_model_path="models",
        ms2_tolerance=0.05,
        calibration_set_size=0.15,
        skip_deeplc_retrain=True,
        processes=2,
        id_decoy_pattern="^DECOY_",
        lower_score_is_better=True,
        log_level="INFO",
        spectrum_id_pattern="(.*)",
        psm_id_pattern="(.*)",
    )
    annotator.build_idxml_data(idxml_file, mzml_file)
    annotator.annotate()

    output_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet_rescored.idXML"
    )

    annotator.write_idxml_file(output_file)


def test_idxmlreader_wrong_model():

    idxml_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML"
    )

    mzml_file = (
        TESTS_DIR / "test_data" / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML"
    )

    annotator = Annotator(
        feature_generators="ms2pip,deeplc",
        only_features="DeepLC:RtDiff,DeepLC:PredictedRetentionTimeBest,Ms2pip:DotProd",
        ms2pip_model="CID-TMT",
        ms2pip_model_path="models",
        ms2_tolerance=0.05,
        calibration_set_size=0.15,
        skip_deeplc_retrain=True,
        processes=2,
        id_decoy_pattern="^DECOY_",
        lower_score_is_better=True,
        log_level="INFO",
        spectrum_id_pattern="(.*)",
        psm_id_pattern="(.*)",
        find_best_ms2pip_model=True,
    )
    annotator.build_idxml_data(idxml_file, mzml_file)
    annotator.annotate()

    output_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet_rescored.idXML"
    )

    annotator.write_idxml_file(output_file)


@pytest.mark.skip(reason="This is for local test in big datasets, kipping for now")
def test_idxmlreader_failing_help():
    idxml_file = (
        TESTS_DIR
        / "test_data"
        / "dae1cb16fb57893b94bfcb731b2bf7"
        / "Instrument1_sample14_S1R10_042116_Fr12_msgf.idXML"
    )

    mzml_file = (
        TESTS_DIR
        / "test_data"
        / "dae1cb16fb57893b94bfcb731b2bf7"
        / "Instrument1_sample14_S1R10_042116_Fr12.mzML"
    )

    idexml_reader = IdXMLRescoringReader(idexml_filename=idxml_file)
    peptide_list = idexml_reader._build_psm_index()

    psm_count = OpenMSHelper.get_psm_count(idexml_reader.oms_peptides)

    logging.info("Loaded %s PSMs from %s", psm_count)
    assert len(peptide_list) == 66770

    openms_helper = OpenMSHelper()
    decoys, targets = openms_helper.count_decoys_targets(idexml_reader.oms_peptides)
    logging.info(
        "Loaded %s PSMs from %s, %s decoys and %s targets", psm_count, idxml_file, decoys, targets
    )
    assert decoys == 25171
    assert targets == 41599

    idexml_reader._build_spectrum_lookup(mzml_file)
    missing_count = idexml_reader.validate_psm_spectrum_references()

    assert missing_count == 0


def test_sage_feature_file():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "sage2feature",
            "--idxml",
            f"{TESTS_DIR}/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_sage_ms2rescore.idXML",
            "--output_file",
            f"{TESTS_DIR}/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_sage_ms2rescore_feat_gen.idXML",
            "--feat_file",
            f"{TESTS_DIR}/test_data/tmt_erwinia_1ulsike_top10hcd_isol2_45stepped_60min_01_sage_ms2rescore.idxml.feature_names.tsv",
        ],
    )

    assert result.exit_code == 0


def test_spectrum2feature_file():
    idxml_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_sage_ms2rescore.idXML"
    )
    mzml_file = (
        TESTS_DIR / "test_data" / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML"
    )
    output_file = (
        TESTS_DIR
        / "test_data"
        / "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_sage_ms2rescore_snr.idXML"
    )
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "spectrum2feature",
            "--mzml",
            mzml_file,
            "--idxml",
            idxml_file,
            "--output",
            output_file,
        ],
    )

    assert result.exit_code == 0


# test for the add_sage_feature command in cli
def test_add_sage_feature_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["sage2feature", "--help"])

    assert result.exit_code == 0


# test for the ms2rescore


@pytest.mark.skip(reason="This is for local test in big datasets, kipping for now")
def test_local_file():
    runner = CliRunner()
    local_folder = TESTS_DIR / "test_data" / "dae1cb16fb57893b94bfcb731b2bf7"
    result = runner.invoke(
        cli,
        [
            "msrescore2feature",
            "--idxml",
            "{}/UPS1_50amol_R1_comet.idXML".format(local_folder),
            "--mzml",
            "{}/UPS1_50amol_R1_converted.mzML".format(local_folder),
            "--processes",
            "2",
            "--ms2pip_model",
            "CID",
            "--feature_generators",
            "ms2pip,deeplc",
            "--id_decoy_pattern",
            "^DECOY_",
            "--output",
            "{}/UPS1_50amol_R1_rescored.idXML".format(local_folder),
            "--ms2_tolerance",
            "0.05",
        ],
    )
    assert result.exit_code == 0
