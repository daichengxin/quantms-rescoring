import click


from quantmsrescore import __version__
from quantmsrescore.ms2rescore import ms2rescore

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.version_option(
    version=__version__, package_name="quantmsrescore", message="%(package)s %(version)s"
)
@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(ms2rescore)


def main():
    try:
        cli()
    except SystemExit as e:
        if e.code != 0:
            raise

if __name__ == "__main__":
    main()
