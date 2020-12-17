"""Console script for ravenpy."""
import json
import sys

import click

from ravenpy.models.importers import RoutingProductGridWeightImporter


@click.group()
def cli():
    pass


@cli.command(help="Generate grid weights in various formats.")
@click.argument("shp-file-path", type=click.Path(exists=True))
@click.argument("nc-file-path", type=click.Path(exists=True))
@click.option(
    "-d",
    "--dim-names",
    nargs=2,
    type=str,
    default=RoutingProductGridWeightImporter.DIM_NAMES,
    show_default=True,
)
@click.option(
    "-v",
    "--var-names",
    nargs=2,
    type=str,
    default=RoutingProductGridWeightImporter.VAR_NAMES,
    show_default=True,
)
@click.option(
    "-o",
    "--output",
    type=click.Choice(["raven", "json", "text"]),
    default="raven",
    show_default=True,
)
def generate_grid_weights(shp_file_path, nc_file_path, dim_names, var_names, output):
    importer = RoutingProductGridWeightImporter(
        shp_file_path, nc_file_path, dim_names, var_names
    )
    cmd = importer.extract()
    if output == "raven":
        print(cmd.to_rv())
    elif output == "text":
        print("\n".join(" ".join(map(str, d)) for d in cmd.data))
    elif output == "json":
        print(
            json.dumps(
                {
                    "NumberHRUs": cmd.number_hrus,
                    "NumberGridCells": cmd.number_grid_cells,
                    "Weights": cmd.data,
                },
                indent=4,
            )
        )


if __name__ == "__main__":
    sys.exit(cli())  # pragma: no cover
