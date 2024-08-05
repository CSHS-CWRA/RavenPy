import click

from ravenpy import __version__

from .aggregate_forcings_to_hrus import aggregate_forcings_to_hrus
from .collect_subbasins_upstream_of_gauge import collect_subbasins_upstream_of_gauge
from .generate_grid_weights import generate_grid_weights
from .generate_hrus_from_routing_product import generate_hrus_from_routing_product


@click.group()
@click.version_option(__version__)
def cli():
    pass


def main():
    cli.add_command(generate_grid_weights)
    cli.add_command(aggregate_forcings_to_hrus)
    cli.add_command(collect_subbasins_upstream_of_gauge)
    cli.add_command(generate_hrus_from_routing_product)
    cli()
