from ravenpy.models.importers import (
    RoutingProductGridWeightImporter,  
)
import netCDF4 as nc4
import numpy as np


DIM_NAMES = ("lon", "lat")
VAR_NAMES = ("lon", "lat")
ROUTING_ID_FIELD = "HRU_ID"
NETCDF_INPUT_FIELD = "NetCDF_col"
AREA_ERROR_THRESHOLD = 0.05
VARIABLES_TO_AGGREGATE = ["Streaminputs"]

def aggregate_forcings_to_HRUs( input_file, routing_file, output_file,
                                    dim_names=DIM_NAMES,
                                    var_names=VAR_NAMES,
                                    routing_id_field=ROUTING_ID_FIELD,
                                    netcdf_input_field=NETCDF_INPUT_FIELD,
                                    gauge_ids=None,
                                    sub_ids=None,
                                    area_error_threshold=AREA_ERROR_THRESHOLD,
                                    variables_to_aggregate=VARIABLES_TO_AGGREGATE,
                                    ):

    """
    arguments:
         all same as in RoutingProductGridWeightImporter
         variables_to_aggregate ... list of strings of netcdf variables that will be aggregated
                                    variable needs to be 3D (lat, lon, time) - dimension order does not matter
    """

    
    importer = RoutingProductGridWeightImporter(input_file, routing_file,
                                                    dim_names=dim_names, var_names=var_names,
                                                    routing_id_field=routing_id_field,
                                                    netcdf_input_field=netcdf_input_field,
                                                    gauge_ids=gauge_ids,
                                                    sub_ids=sub_ids,
                                                    area_error_threshold=area_error_threshold)
    gws = importer.extract()

    nHRU         = gws.number_hrus
    nCells       = gws.number_grid_cells
    weights_data = gws.data


    # read NetCDF
    nc_in = nc4.Dataset(input_file,'r')

    # length of dimensions
    nlon  = nc_in.dimensions[dim_names[0]].size
    nlat  = nc_in.dimensions[dim_names[1]].size
    ntime = nc_in.dimensions['time'].size

    # check that this matches the nCells from grid-weights
    # (basically making sure these are the same NetCDF files; not sure how that would happen but just to be on the safe side)
    if nlon*nlat != nCells:
        raise ValueError(
                "Number of cells used to derive grid-weights does not match NetCDF provided for aggregating variables"
            )

    # convert cell_id (from weights) into (lat_id and lon_id)
    # cell_id = ilat * nlon + ilon
    # ---> ilon = cell_id %  nlon
    # ---> ilat = cell_id // nlon
    weights_data_lon_lat_ids = []
    for iweights_data in weights_data:

        cell_id = iweights_data[1]
        ilon = cell_id %  nlon
        ilat = cell_id // nlon
        weights_data_lon_lat_ids.append( tuple( [ iweights_data[0],ilon,ilat,iweights_data[2] ] ) )
    weights_data_lon_lat_ids = np.array(weights_data_lon_lat_ids)   # because I dont know how to do most operations elegantly with lists of tuples

    # create new NetCDF that will contain aggregated data of listed variables
    nc_out = nc4.Dataset(output_file,'w')
    nc_dim_time = nc_out.createDimension('time',ntime)
    nc_dim_hrus = nc_out.createDimension('nHRU',nHRU)

    # copy all global attributes over
    nc_out.setncatts(nc_in.__dict__)

    # create all variables in output NC (incl. time) and copy over all attributes
    for name, variable in nc_in.variables.items():
        if name in variables_to_aggregate+['time']:

            if name != 'time':
                dims = ['time','nHRU']
            else:
                dims = ['time']
            x = nc_out.createVariable(name, variable.datatype, dims, zlib=True)

            # copy variable attributes all at once via dictionary
            nc_out[name].setncatts(nc_in[name].__dict__)

            # copy over time variable data
            if name == 'time':
                nc_out[name][:] = nc_in[name][:]


    for variable_to_aggregate in variables_to_aggregate:

        # read 3D variable
        input_var  = nc_in.variables[variable_to_aggregate]

        # read 2D variable (to be able to write to it later)
        output_var  = nc_out.variables[variable_to_aggregate]

        # is variable 3D ?
        ndims = len(input_var.dimensions)
        if ndims != 3:
            raise ValueError(
                    "Input variable to aggregate needs to be 3D"
                )
        
        # what is the order of dimensions?
        idx_lon_dim  = input_var.dimensions.index( dim_names[0] )
        idx_lat_dim  = input_var.dimensions.index( dim_names[1] )
        idx_time_dim = input_var.dimensions.index('time')

        # list of HRU IDs
        hrus = np.unique( weights_data_lon_lat_ids[:,0] )

        if len(hrus) != nHRU:
            # should really never happen
            raise ValueError(
                "Number of weights found in grid weights list is not matching the number indicated there by nHRUs"
                )

        # do actual aggregation
        agg_var = np.zeros( [ntime,nHRU] )
        new_weights = []
        for ihru,hru in enumerate(hrus):

            hru=int(hru)

            # filter all weights for current HRU
            idx = np.where( weights_data_lon_lat_ids[:,0] == hru )[0]

            for ii in idx:
                # bring idx for input_var in appropriate order
                idx_input = [slice(0,ntime,1), slice(0,ntime,1), slice(0,ntime,1)]
                idx_input[idx_lon_dim] = int(weights_data_lon_lat_ids[ii,1])
                idx_input[idx_lat_dim] = int(weights_data_lon_lat_ids[ii,2])
                agg_var[:,ihru] +=  input_var[idx_input] * weights_data_lon_lat_ids[ii,3]

            # create new weights: now each HRU is exactly one "grid-cell"
            new_weights.append( tuple( [hru,ihru,1.0] ) )

        # write 2D variable
        output_var[:] = agg_var
            
    nc_out.close()

    # the return should actually look exactly like the return of "RoutingProductGridWeightImporter"
    # not sure how to do that

    return new_weights





# input_file = "/Users/j6mai/Documents/GitHub/GridWeightsGenerator/example/input_VIC/VIC_streaminputs.nc"	
# routing_file = "/Users/j6mai/Documents/GitHub/GridWeightsGenerator/example/maps/finalcat_hru_info.zip"
# output_file = "/Users/j6mai/Documents/GitHub/GridWeightsGenerator/example/input_VIC/VIC_streaminputs_agg.nc"	
# variables_to_aggregate=VARIABLES_TO_AGGREGATE
# dim_names=DIM_NAMES
# var_names=VAR_NAMES
# routing_id_field=ROUTING_ID_FIELD
# netcdf_input_field=NETCDF_INPUT_FIELD
# gauge_ids=None
# sub_ids=None
# area_error_threshold=AREA_ERROR_THRESHOLD




