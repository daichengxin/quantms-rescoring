# Handle environment variable before any imports
import os

# Temporarily unset OPENMS_DATA_PATH to prevent warning
openms_data_path = os.environ.get("OPENMS_DATA_PATH")
if openms_data_path:
    del os.environ["OPENMS_DATA_PATH"]

# Import modules
import click
from quantmsrescore import __version__
from quantmsrescore.ms2rescore import msrescore2feature
from quantmsrescore.sage_feature import add_sage_feature
from quantmsrescore.snr import spectrum2feature
from quantmsrescore.psm_clean import psm_feature_clean

# Restore the environment variable
if openms_data_path:
    os.environ["OPENMS_DATA_PATH"] = openms_data_path

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.version_option(
    version=__version__, package_name="quantmsrescore", message="%(package)s %(version)s"
)
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(msrescore2feature)
cli.add_command(add_sage_feature)
cli.add_command(spectrum2feature)
cli.add_command(psm_feature_clean)


def main():
    try:
        cli()
    except SystemExit as e:
        if e.code != 0:
            raise


if __name__ == "__main__":
    main()
