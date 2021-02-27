import click

from ravenpy import __version__

from .generate_grid_weights import generate_grid_weights


@click.group()
@click.version_option(__version__)
def cli():
    pass


def main():
    cli.add_command(generate_grid_weights)
    cli()
