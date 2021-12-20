from pathlib import Path

import click

from ravenpy.extractors.routing_product import RoutingProductGridWeightExtractor


@click.command()
@click.argument("input-file", type=click.Path(exists=True))
@click.argument("routing-file", type=click.Path(exists=True))
@click.option(
    "-d",
    "--dim-names",
    nargs=2,
    type=click.Tuple([str, str]),
    default=RoutingProductGridWeightExtractor.DIM_NAMES,
    show_default=True,
    help="Ordered dimension names of longitude (x) and latitude (y) in the NetCDF INPUT_FILE.",
)
@click.option(
    "-v",
    "--var-names",
    nargs=2,
    type=click.Tuple([str, str]),
    default=RoutingProductGridWeightExtractor.VAR_NAMES,
    show_default=True,
    help="Variable name of 1D or 2D longitude and latitude variables in the NetCDF INPUT_FILE (in this order).",
)
@click.option(
    "-c",  # legacy script short option
    "--routing-id-field",
    type=str,
    default=RoutingProductGridWeightExtractor.ROUTING_ID_FIELD,
    show_default=True,
    help="Name of column in routing information shapefile (ROUTING_FILE) containing a unique key for each dataset.",
)
@click.option(
    "-f",  # legacy script short option
    "--netcdf-input-field",
    type=str,
    default=RoutingProductGridWeightExtractor.NETCDF_INPUT_FIELD,
    show_default=True,
    help="Attribute name in INPUT_FILE shapefile that defines the index of the shape in NetCDF model output file (numbering needs to be [0 ... N-1]).",
)
@click.option(
    "-g",
    "--gauge-id",
    "gauge_ids",  # rename arg to put emphasis on the fact that it's a list
    multiple=True,
    type=str,
    show_default=True,
    help="Streamflow gauge IDs of interest (corresponds to 'Gauge_ID' in the ROUTING_FILE shapefile).",
)
@click.option(
    "-s",
    "--sub-id",
    "sub_ids",  # rename arg to put emphasis on the fact that it's a list
    multiple=True,
    type=int,
    show_default=True,
    help="IDs of subbasins of interest (containing usually a gauge station, corresponds to 'SubId' in the ROUTING_FILE shapefile).",
)
@click.option(
    "-e",  # legacy script short option
    "--area-error-threshold",
    type=float,
    default=RoutingProductGridWeightExtractor.AREA_ERROR_THRESHOLD,
    show_default=True,
    help="Threshold (as fraction) of allowed mismatch in areas between subbasins from routing information (ROUTING_FILE) and overlay with grid-cells or subbasins (INPUT_FILE). If error is smaller than this threshold the weights will be adjusted such that they sum up to exactly 1. Raven will exit gracefully in case weights do not sum up to at least 0.95.",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    help="Text field that will contain the results as a single :GridWeights Raven command containing the weights.",
)
def generate_grid_weights(
    input_file,
    routing_file,
    dim_names,
    var_names,
    routing_id_field,
    netcdf_input_field,
    gauge_ids,
    sub_ids,
    area_error_threshold,
    output,
):
    """
    Generate grid weights in various formats.

    INPUT_FILE: File containing model discretization. Can be either:

    (A) NetCDF file containing at least 1D or 2D latitudes and 1D or 2D longitudes where this
    grid needs to be representative of model outputs that are then required to be routed. The names of the dimensions
    and the variables holding the lat/lon information should be specified with options --dim-names (-d) and --var-names (-v).

    (B) Shapefile (either a .shp or a .zip) that contains shapes of subbasins and one attribute in this
    shapefile that is defining its index in the NetCDF model output file (numbering needs to be [0 ... N-1]). The name of
    this attribute should be specified via option --netcdf-input-field (-f).

    ROUTING_FILE: Shapefile (either a .shp or a .zip) that contains all information of the routing toolbox for the catchment
    of interest (and maybe some more catchments). This file should contain an attribute with a unique identifier for the HRUs:
    the default is "HRU_ID", but it can be set with --routing-id-field (-c).

    The script will output the results as RVT file with a single :GridWeights command block containing the weights (that the user
    is then free to embed or reference, in her own config context).
    """
    # NOTE: This is in order to make sphinx-click happy. Magic. Do not touch.
    from ravenpy.extractors import RoutingProductGridWeightExtractor

    extractor = RoutingProductGridWeightExtractor(
        input_file,
        routing_file,
        dim_names,
        var_names,
        routing_id_field,
        netcdf_input_field,
        gauge_ids,
        sub_ids,
        area_error_threshold,
    )
    gw_cmd = extractor.extract()

    if not output:
        input_file_path = Path(input_file)
        output_file_path = (
            input_file_path.parent / f"{input_file_path.stem}_weights.rvt"
        )
    else:
        output_file_path = Path(output)

    output_file_path.write_text(gw_cmd.to_rv() + "\n")

    click.echo(f"Created {output_file_path}")
