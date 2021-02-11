import numpy as np
import netCDF4 as nc4
import pytest

from ravenpy.models.aggregate_forcings_to_HRUs import (
    aggregate_forcings_to_HRUs,
    convert_cell_id_to_ilat_ilon
    )
from ravenpy.utilities.testdata import get_local_testdata

from ravenpy.models.importers import (
    RoutingProductGridWeightImporter,
    GridWeightsCommand
    )


# run test:
#   rpy
#   cd /Users/j6mai/Documents/GitHub/RavenPy/
#   pytest tests/test_aggregate_forcings_to_HRUs.py


def test_convert_cell_id_to_ilat_ilon():

    input_file = get_local_testdata("raven-routing-sample/VIC_streaminputs.nc")
    routing_file = get_local_testdata("raven-routing-sample/finalcat_hru_info.zip")

    DIM_NAMES = ("lon", "lat")
    VAR_NAMES = ("lon", "lat")
    ROUTING_ID_FIELD = "HRU_ID"
    NETCDF_INPUT_FIELD = "NetCDF_col"
    AREA_ERROR_THRESHOLD = 0.05

    importer = RoutingProductGridWeightImporter(input_file, routing_file, 
                                    dim_names=DIM_NAMES,
                                    var_names=VAR_NAMES,
                                    routing_id_field=ROUTING_ID_FIELD,
                                    netcdf_input_field=NETCDF_INPUT_FIELD,
                                    gauge_ids=None,
                                    sub_ids=None,
                                    area_error_threshold=AREA_ERROR_THRESHOLD)
    gws = importer.extract()
    old_ids = gws.data

    # read dimensions from NetCDF
    nc_in = nc4.Dataset(input_file,'r')
    nlon  = nc_in.dimensions[DIM_NAMES[0]].size
    nlat  = nc_in.dimensions[DIM_NAMES[1]].size

    gws_new = convert_cell_id_to_ilat_ilon(gws, nlon, nlat)
    new_ids = gws_new.data

    assert new_ids[4][1] == 4.0                                # lon_id   5th grid weight
    assert new_ids[4][2] == 5.0                                # lat_id   5th grid weight
    assert new_ids[4][2]*nlon + new_ids[4][1] == old_ids[4][1] # derived cell_id --> ilat * nlon + ilon
    assert old_ids[4][1] == 54                                 # cell_id  5th grid weight
    

    


def test_aggregate_forcings_to_HRUs():

    input_file = get_local_testdata("raven-routing-sample/VIC_streaminputs.nc")
    routing_file = get_local_testdata("raven-routing-sample/finalcat_hru_info.zip")
    output_file = "/tmp/VIC_streaminputs_aggregated.nc"   # this should be something better

    DIM_NAMES = ("lon", "lat")
    VAR_NAMES = ("lon", "lat")
    ROUTING_ID_FIELD = "HRU_ID"
    NETCDF_INPUT_FIELD = "NetCDF_col"
    AREA_ERROR_THRESHOLD = 0.05
    VARIABLES_TO_AGGREGATE = ["Streaminputs"]

    gws = aggregate_forcings_to_HRUs( input_file, routing_file, output_file,
                                    dim_names=DIM_NAMES,
                                    var_names=VAR_NAMES,
                                    routing_id_field=ROUTING_ID_FIELD,
                                    netcdf_input_field=NETCDF_INPUT_FIELD,
                                    gauge_ids=None,
                                    sub_ids=None,
                                    area_error_threshold=AREA_ERROR_THRESHOLD,
                                    variables_to_aggregate=VARIABLES_TO_AGGREGATE )

    new_weights = gws.data

    # check new weights
    assert new_weights[0][0] == 1            # These are the HRU-IDs
    
    assert new_weights[0][1] == 0            # These need to be exactly [0,1,2,3,...,nHRU]
    assert new_weights[1][1] == 1            # These need to be exactly [0,1,2,3,...,nHRU]
    assert new_weights[2][1] == 2            # These need to be exactly [0,1,2,3,...,nHRU]
    assert new_weights[3][1] == 3            # These need to be exactly [0,1,2,3,...,nHRU]
    
    assert new_weights[0][2] == 1.0          # All new_weights[:][2] need to be 1.0
    assert new_weights[1][2] == 1.0          # All new_weights[:][2] need to be 1.0
    assert new_weights[2][2] == 1.0          # All new_weights[:][2] need to be 1.0
    assert new_weights[3][2] == 1.0          # All new_weights[:][2] need to be 1.0


    # check aggregated NetCDF file 
    nc_in = nc4.Dataset(output_file,'r')
    val = nc_in.variables[VARIABLES_TO_AGGREGATE[0]][:]
    nc_in.close()
    assert abs(val[    0, 0] - 0.017309) < 1e-04
    assert abs(val[16071, 0] - 0.569977) < 1e-04
    assert abs(val[    0,50] - 0.010276) < 1e-04
    assert abs(val[16071,50] - 0.516639) < 1e-04
