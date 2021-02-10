from ravenpy.models.importers import (
    RoutingProductGridWeightImporter,  
)
import netCDF4 as nc4


DIM_NAMES = ("lon", "lat")
VAR_NAMES = ("lon", "lat")
ROUTING_ID_FIELD = "HRU_ID"
NETCDF_INPUT_FIELD = "NetCDF_col"
AREA_ERROR_THRESHOLD = 0.05
VARIABLES_TO_AGGREGATE = ["Streaminputs"]

def aggregate_forcings_to_HRUs( input_file, routing_file,
                                    dim_names=DIM_NAMES,
                                    var_names=VAR_NAMES,
                                    routing_id_field=ROUTING_ID_FIELD,
                                    netcdf_input_field=NETCDF_INPUT_FIELD,
                                    gauge_ids=None,
                                    sub_ids=None,
                                    area_error_threshold=AREA_ERROR_THRESHOLD,
                                    variables_to_aggregate=VARIABLES_TO_AGGREGATE ):

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
    
    # cell_id = ilat * self._nlon + ilon
    print(weights_data[20])
    
    
    # read NetCDF
    input_data = nc4.Dataset(input_file)
    
    for variable_to_aggregate in variables_to_aggregate:
        input_var  = input_data.variables[variable_to_aggregate]

    return nHRU    # dummy return for now :)



# from ravenpy.utilities.testdata import get_local_testdata

# input_file = get_local_testdata("raven-routing-sample/VIC_streaminputs.nc")
# routing_file = get_local_testdata("raven-routing-sample/finalcat_hru_info.zip")

# aggregate_forcings_to_HRUs( input_file, routing_file,
#                                     dim_names=DIM_NAMES,
#                                     var_names=VAR_NAMES,
#                                     routing_id_field=ROUTING_ID_FIELD,
#                                     netcdf_input_field=NETCDF_INPUT_FIELD,
#                                     gauge_ids=None,
#                                     sub_ids=None,
#                                     area_error_threshold=AREA_ERROR_THRESHOLD,
#                                     variables_to_aggregate=VARIABLES_TO_AGGREGATE )
