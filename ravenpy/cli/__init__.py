import click

from ravenpy import __version__

from .aggregate_forcings_to_hrus import aggregate_forcings_to_hrus
from .generate_grid_weights import generate_grid_weights


@click.group()
@click.version_option(__version__)
def cli():
    pass


def main():
    cli.add_command(generate_grid_weights)
    cli.add_command(aggregate_forcings_to_hrus)
    cli()
