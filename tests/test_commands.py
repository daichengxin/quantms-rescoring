from click.testing import CliRunner

from quantmsrescore.rescoring import cli
from pathlib import Path
import logging

TESTS_DIR = Path(__file__).parent

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# test for the rescoring command in cli
def test_ms2rescore_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["ms2rescore", "--help"])

    assert result.exit_code == 0


def test_ms2rescore():
    runner = CliRunner()
    result = runner.invoke(
        cli,
        [
            "ms2rescore",
            "--psm_file",
            "{}/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML".format(TESTS_DIR),
            "--spectrum_path",
            "{}/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML".format(TESTS_DIR),
            "--processes",
            "2",
            "--ms2pip_model",
            "HCD2021",
            "--feature_generators",
            "'ms2pip,deeplc'",
            "--id_decoy_pattern",
            "^rev",
            "--test_fdr",
            "0.05",
        ],
    )
    assert result.exit_code == 0
