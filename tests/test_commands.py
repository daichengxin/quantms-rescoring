from click.testing import CliRunner

from quantmsrescore.rescoring import cli


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
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01_comet.idXML",
            "--spectrum_path",
            "tests/test_data/TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML",
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
