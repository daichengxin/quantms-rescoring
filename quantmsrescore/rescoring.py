# Handle environment variables and warnings before any imports
import os
import sys
import warnings

# Suppress all warnings related to OPENMS_DATA_PATH
warnings.filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*", category=UserWarning)
warnings.filterwarnings("ignore", message=".*OPENMS_DATA_PATH.*")
warnings.filterwarnings("ignore", category=UserWarning, module="pyopenms")

# Store the original stdout to restore it later
original_stdout = sys.stdout
# Temporarily redirect stdout to suppress any output during import
sys.stdout = open(os.devnull, 'w')

# Try to handle the environment variable directly
openms_data_path = os.environ.get("OPENMS_DATA_PATH")
if openms_data_path:
    # Temporarily unset the variable
    del os.environ["OPENMS_DATA_PATH"]

# After all our preparations, we can now safely import other modules
import click

from quantmsrescore import __version__
from quantmsrescore.ms2rescore import msrescore2feature
from quantmsrescore.sage_feature import add_sage_feature
from quantmsrescore.snr import spectrum2feature

# Restore the environment variable if we unset it
if openms_data_path:
    os.environ["OPENMS_DATA_PATH"] = openms_data_path

# Restore stdout
sys.stdout.close()
sys.stdout = original_stdout

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


def main():
    try:
        cli()
    except SystemExit as e:
        if e.code != 0:
            raise


if __name__ == "__main__":
    main()
