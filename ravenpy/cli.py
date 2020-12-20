"""Console script for ravenpy."""
import json
import re
import sys

import click

from ravenpy.models.importers import RoutingProductGridWeightImporter


@click.group()
def cli():
    pass


@cli.command()
@click.argument("shp-file-path", type=click.Path(exists=True))
@click.argument("nc-file-path", type=click.Path(exists=True))
@click.option(
    "-d",
    "--dim-names",
    nargs=2,
    type=click.Tuple([str, str]),
    default=RoutingProductGridWeightImporter.DIM_NAMES,
    show_default=True,
    help="Ordered dimension names of longitude (x) and latitude (y) in the NetCDF file.",
)
@click.option(
    "-v",
    "--var-names",
    nargs=2,
    type=click.Tuple([str, str]),
    default=RoutingProductGridWeightImporter.VAR_NAMES,
    show_default=True,
    help="Ordered variable names of longitude (x) and latitude (y) in the NetCDF file.",
)
@click.option(
    "--hru-id-field",
    type=str,
    default=RoutingProductGridWeightImporter.HRU_ID_FIELD,
    show_default=True,
)
@click.option(
    "-g",
    "--gauge-id",
    multiple=True,
    type=str,
    show_default=True,
    help="Basins of interest (corresponds to 'Gauge-ID' in the shapefile).",
)
@click.option(
    "-s",
    "--sub-id",
    multiple=True,
    type=str,  # to support alternate "-s 123,456" legacy syntax
    show_default=True,
    help="Sub IDs of most downstream subbasins (containing usually a gauge station, corresponds to 'SubId' in the shapefile).",
)
@click.option(
    "-o",
    "--output",
    type=click.Choice(["raven", "json", "text"]),
    default="raven",
    show_default=True,
)
def generate_grid_weights(
    shp_file_path,
    nc_file_path,
    dim_names,
    var_names,
    hru_id_field,
    gauge_id,
    sub_id,
    output,
):
    """
    Generate grid weights in various formats.

    SHP_FILE_PATH: Shapefile that contains all information of the routing toolbox for the catchment of interest (and maybe some more catchments).

    NC_FILE_PATH: NetCDF file containing at least 2D latitudes and 2D longitudes. Grid needs to be representative of model outputs that are then required to be routed.
    """

    # Although Click does not support it, we want to support legacy syntax (e.g. "-s 123,456") for both gauge and sub IDs.
    # The Click way for `multiple` options would be: "-s 123 -s 456"

    if len(gauge_id) == 1:
        gauge_ids = list(filter(None, re.split("\W+", gauge_id[0])))
    else:
        gauge_ids = gauge_id

    if len(sub_id) == 1:
        sub_ids = list(map(int, filter(None, re.split("\W+", sub_id[0]))))
    else:
        sub_ids = list(map(int, sub_id))

    importer = RoutingProductGridWeightImporter(
        shp_file_path,
        nc_file_path,
        dim_names,
        var_names,
        hru_id_field,
        gauge_ids,
        sub_ids,
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
