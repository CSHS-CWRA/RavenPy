"""Console script for ravenpy."""
import json
import re
import sys

import click

from ravenpy import __version__
from ravenpy.models.importers import RoutingProductGridWeightImporter


@click.group()
@click.version_option(__version__)
def cli():
    pass


@cli.command()
@click.argument("input-file-path", type=click.Path(exists=True))
@click.argument("routing-file-path", type=click.Path(exists=True))
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
    help="(A) Variable name of 2D longitude and latitude variables in NetCDF (in this order). Example: 'lon,lat'. Or (B) Attribute name in input_file shapefile (option -i) that defines the index of the shape in NetCDF model output file (numbering needs to be [0 ... N-1]).",
)
@click.option(
    "-c",  # legacy script short option
    "--routing-id-field",
    type=str,
    default=RoutingProductGridWeightImporter.ROUTING_ID_FIELD,
    show_default=True,
    help="Name of column in routing information shapefile containing unique key for each dataset.",
)
@click.option(
    "-f",  # legacy script short option
    "--netcdf-input-field",
    type=str,
    default=RoutingProductGridWeightImporter.NETCDF_INPUT_FIELD,
    show_default=True,
    help="Attribute name in input_file shapefile that defines the index of the shape in NetCDF model output file (numbering needs to be [0 ... N-1]).",
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
    "-e",  # legacy script short option
    "--area-error-threshold",
    type=float,
    default=RoutingProductGridWeightImporter.AREA_ERROR_THRESHOLD,
    show_default=True,
    help="Threshold (as fraction) of allowed mismatch in areas between subbasins from routing information (-r) and overlay with grid-cells or subbasins (-i). If error is smaller than this threshold the weights will be adjusted such that they sum up to exactly 1. Raven will exit gracefully in case weights do not sum up to at least 0.95.",
)
@click.option(
    "-o",
    "--output",
    type=click.Choice(["raven", "json", "text"]),
    default="raven",
    show_default=True,
)
def generate_grid_weights(
    input_file_path,
    routing_file_path,
    dim_names,
    var_names,
    routing_id_field,
    netcdf_input_field,
    gauge_id,
    sub_id,
    area_error_threshold,
    output,
):
    """
    Generate grid weights in various formats.

    INPUT_FILE_PATH: Either (A) Example NetCDF file containing at least 1D or 2D latitudes and 1D or 2D longitudes where grid
    needs to be representative of model outputs that are then required to be routed. Or (B) a shapefile (either a .shp or a .zip)
    that contains shapes of subbasins and one attribute that is defining its index in the NetCDF model output file (numbering
    needs to be [0 ... N-1]).

    ROUTING_FILE_PATH: Shapefile (either a .shp or a .zip) that contains all information of the routing toolbox for the catchment
    of interest (and maybe some more catchments).

    The script outputs the results (in the chosen format) on STDOUT, so they can be redirected to a file using the `>` shell operator.
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
        input_file_path,
        routing_file_path,
        dim_names,
        var_names,
        routing_id_field,
        netcdf_input_field,
        gauge_ids,
        sub_ids,
        area_error_threshold,
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
