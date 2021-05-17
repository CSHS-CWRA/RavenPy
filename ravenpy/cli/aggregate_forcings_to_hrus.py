from pathlib import Path
from typing import List, Tuple

import click

from ravenpy.config.commands import GridWeightsCommand
from ravenpy.extractors.routing_product import RoutingProductGridWeightExtractor


@click.command()
@click.argument("input-nc-file", type=click.Path(exists=True))
@click.argument("input-weight-file", type=click.Path(exists=True))
@click.option(
    "-d",
    "--dim-names",
    nargs=2,
    type=click.Tuple([str, str]),
    default=RoutingProductGridWeightExtractor.DIM_NAMES,
    show_default=True,
    help="Ordered dimension names of longitude (x) and latitude (y) in the NetCDF INPUT_NC_FILE.",
)
@click.option(
    "-v",
    "--var-to-aggregate",
    "variables_to_aggregate",  # rename arg to put emphasis on the fact that it's a list
    required=True,
    multiple=True,
    type=str,
    show_default=True,
    help="Variables to aggregate in INPUT_NC_FILE (at least one).",
)
@click.option("--output-nc-file", type=click.Path(), help="")
@click.option("--output-weight-file", type=click.Path(), help="")
def aggregate_forcings_to_hrus(
    input_nc_file,
    input_weight_file,
    output_nc_file,
    output_weight_file,
    dim_names,
    variables_to_aggregate,
):
    """
    Aggregates NetCDF files containing 3-dimensional forcing variables like precipitation and temperature
    over (x,y,time) into 2-dimensional forcings for each of the n HRUs of a specific basin over (n,time).
    The 3-dimensional NetCDF files are usually used in ``:GriddedForcing`` commands in Raven while the 2-dimensional
    ones can be used in ``:StationForcing`` commands. The NetCDF files generated with this function will only
    contain the forcings required to simulate an individual basin and hence file sizes are smaller and Raven
    runtimes can decrease drastically under certain conditions.

    INPUT_NC_FILE: NetCDF file containing 3-dimensional variables that will be aggregated. Either all variables
    will be aggregated or only a subset specified using --var-to-aggregate (e.g., [precip,temp]). The name of
    the spatial dimensions of the NetCDF are assumed to be (lon_dim, lat_dim). Otherwise they will need to be
    specified using --dim-names. The order of the three dimensions for each variable does not matter; the
    function will arrange them as required.

    INPUT_WEIGHT_FILE: A text file containing the grid weights derived using the script "generate-grid-weights"
    for the basin forcings are required and the specified NetCDF file. The content of this file must be formatted
    as a valid ``:GridWeights`` Raven command.

    The script outputs two files:

    (1) Aggregated NetCDF file that can be used in a ``:StationForcing`` command in a Raven config.

    (2) A text file (with the same format as INPUT_WEIGHT_FILE) with the updated grid weights, that a ``:StationForcing``
    command will require.
    """
    # NOTE: This is in order to make sphinx-click happy. Magic. Do not touch.
    import netCDF4 as nc4
    import numpy as np

    gws = GridWeightsCommand.parse(Path(input_weight_file).read_text())

    nHRU = gws.number_hrus
    # nCells = gws.number_grid_cells
    weights_data = gws.data

    # read NetCDF
    nc_in = nc4.Dataset(input_nc_file, "r")

    # length of dimensions
    nlon = nc_in.dimensions[dim_names[0]].size
    # nlat = nc_in.dimensions[dim_names[1]].size
    ntime = nc_in.dimensions["time"].size

    # convert weights (cell_id) required by Raven into (lon_id, lat_id)
    # cell_id = ilat * nlon + ilon
    # ---> ilon = cell_id %  nlon
    # ---> ilat = cell_id // nlon
    weights_data_lon_lat_ids = []
    for iweights_data in weights_data:
        cell_id = iweights_data[1]
        ilon = cell_id % nlon
        ilat = cell_id // nlon
        weights_data_lon_lat_ids.append(
            tuple([iweights_data[0], ilon, ilat, iweights_data[2]])
        )
    weights_data_lon_lat_ids = np.asarray(weights_data_lon_lat_ids)  # type: ignore

    # create new NetCDF that will contain aggregated data of listed variables

    if not output_nc_file:
        input_nc_file_path = Path(input_nc_file)
        output_nc_file_path = (
            input_nc_file_path.parent / f"{input_nc_file_path.stem}_aggregated.nc"
        )
    else:
        output_nc_file_path = Path(output_nc_file)

    nc_out = nc4.Dataset(output_nc_file_path, "w")
    _ = nc_out.createDimension("time", ntime)
    _ = nc_out.createDimension("nHRU", nHRU)

    # copy all global attributes over
    nc_out.setncatts(nc_in.__dict__)

    # create all variables in output NC (incl. time) and copy over all attributes
    for name, variable in nc_in.variables.items():
        if name in variables_to_aggregate + ("time",):

            if name != "time":
                dims = ["time", "nHRU"]
            else:
                dims = ["time"]
            _ = nc_out.createVariable(name, variable.datatype, dims, zlib=True)

            # copy variable attributes all at once via dictionary
            nc_out[name].setncatts(nc_in[name].__dict__)

            # copy over time variable data
            if name == "time":
                nc_out[name][:] = nc_in[name][:]

    for variable_to_aggregate in variables_to_aggregate:

        # read 3D variable
        input_var = nc_in.variables[variable_to_aggregate]

        # read 2D variable (to be able to write to it later)
        output_var = nc_out.variables[variable_to_aggregate]

        # is variable 3D ?
        ndims = len(input_var.dimensions)
        if ndims != 3:
            raise ValueError("Input variable to aggregate needs to be 3D")

        # what is the order of dimensions?
        idx_lon_dim = input_var.dimensions.index(dim_names[0])
        idx_lat_dim = input_var.dimensions.index(dim_names[1])
        # idx_time_dim = input_var.dimensions.index("time")

        # read in data for bounding box (is faster than reading then every single cell individually)
        # --> this takes most time for large NetCDFs
        min_lon = int(np.min(weights_data_lon_lat_ids[:, 1]))  # type: ignore
        max_lon = int(np.max(weights_data_lon_lat_ids[:, 1]))  # type: ignore
        min_lat = int(np.min(weights_data_lon_lat_ids[:, 2]))  # type: ignore
        max_lat = int(np.max(weights_data_lon_lat_ids[:, 2]))  # type: ignore

        idx_input = [slice(0, ntime, 1), slice(0, ntime, 1), slice(0, ntime, 1)]
        idx_input[idx_lon_dim] = slice(min_lon, max_lon + 1, 1)
        idx_input[idx_lat_dim] = slice(min_lat, max_lat + 1, 1)

        # print(weights_data_lon_lat_ids[:, 1])
        # print(min_lon, max_lon)
        # print(min_lat, max_lat)
        # print(idx_input)

        input_var_bb = input_var[idx_input]

        # list of HRU IDs
        hrus = np.unique(weights_data_lon_lat_ids[:, 0])  # type: ignore

        if len(hrus) != nHRU:
            # should really never happen
            raise ValueError(
                "Number of weights found in grid weights list is not matching the number indicated there by nHRUs"
            )

        # do actual aggregation
        agg_var = np.zeros([ntime, nHRU])
        new_weights = []
        for ihru, hru in enumerate(hrus):

            hru = int(hru)

            # filter all weights for current HRU
            idx = np.where(weights_data_lon_lat_ids[:, 0] == hru)[0]  # type: ignore

            # create an array that contains weights for each time step
            # --> weights need to be rescaled if NODATA values appear
            # --> transpose such that we can broadcast during rescaling later
            # --> shape = (ncells,ntime) with ncells = len(idx)
            weights_nodata = np.transpose(
                np.tile(weights_data_lon_lat_ids[idx, 3], (ntime, 1))  # type: ignore
            )

            # go through all time steps and zero out weights where grid cell is NODATA
            for iii, ii in enumerate(idx):

                # bring idx for input_var in appropriate order
                idx_input = [slice(0, ntime, 1), slice(0, ntime, 1), slice(0, ntime, 1)]
                idx_input[idx_lon_dim] = int(weights_data_lon_lat_ids[ii, 1]) - min_lon  # type: ignore
                idx_input[idx_lat_dim] = int(weights_data_lon_lat_ids[ii, 2]) - min_lat  # type: ignore

                weights_nodata[iii] = np.where(
                    input_var_bb[idx_input].mask, 0.0, weights_nodata[iii]
                )

            # rescale
            # --> columns where all grid cells are nan might cause a warning:
            #     "RuntimeWarning: invalid value encountered in true_divide"
            old_settings = np.seterr()
            np.seterr(invalid="ignore")
            # "nan_to_num" automatically converts NaN's to 0
            weights_nodata = np.nan_to_num(
                weights_nodata / np.sum(weights_nodata, axis=0)
            )
            np.seterr(**old_settings)  # reset to default

            # derive aggregate
            for iii, ii in enumerate(idx):

                # bring idx for input_var in appropriate order
                idx_input = [slice(0, ntime, 1), slice(0, ntime, 1), slice(0, ntime, 1)]
                idx_input[idx_lon_dim] = int(weights_data_lon_lat_ids[ii, 1]) - min_lon  # type: ignore
                idx_input[idx_lat_dim] = int(weights_data_lon_lat_ids[ii, 2]) - min_lat  # type: ignore

                # this does not work when NODATA values occur
                # agg_var[:, ihru] += (
                #     input_var_bb[idx_input] * weights_data_lon_lat_ids[ii, 3]
                # )

                # this makes sure to handle NODATA cells properly
                agg_var[:, ihru] += input_var_bb[idx_input].data * weights_nodata[iii]

            # if all grid cells have weight 0.0 (i.e., all invalid) set value to NODATA
            idx_t_all_cells_nodata = np.where(
                np.sum(weights_nodata[:, :], axis=0) == 0.0
            )
            if len(idx_t_all_cells_nodata[0]) > 0:
                # only if such cells are found
                # otherwise "input_var.missing_value" might not even exist
                agg_var[idx_t_all_cells_nodata, ihru] = input_var.missing_value

            # create new weights: now each HRU is exactly one "grid-cell"
            new_weights.append((hru, ihru, 1.0))

        # write 2D variable
        output_var[:] = agg_var

    nc_out.close()

    # the return should actually look exactly like the return of "RoutingProductGridWeightImporter"
    # not sure how to do that

    gws_new = GridWeightsCommand(
        number_hrus=nHRU,
        number_grid_cells=nHRU,
        data=tuple(new_weights),
    )

    if not output_weight_file:
        input_weight_file_path = Path(input_weight_file)
        output_weight_file_path = (
            input_weight_file_path.parent
            / f"{input_weight_file_path.stem}_aggregated.rvt"
        )
    else:
        output_weight_file_path = Path(output_weight_file)

    output_weight_file_path.write_text(gws_new.to_rv() + "\n")

    click.echo(f"Created {output_nc_file_path}")
    click.echo(f"Created {output_weight_file_path}")
