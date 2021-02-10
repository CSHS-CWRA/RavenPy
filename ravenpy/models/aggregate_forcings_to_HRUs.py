#!/usr/bin/env python
from __future__ import print_function

# Copyright 2016-2018 Juliane Mai - juliane.mai(at)uwaterloo.ca
#
# License
# This file is part of Juliane Mai's personal code library.
#
# Juliane Mai's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Juliane Mai's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Juliane Mai's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#


from ravenpy.models.importers import (
    RoutingProductGridWeightImporter,
    RoutingProductShapefileImporter,    
)
from ravenpy.models.commands import GriddedForcingCommand
#from ravenpy.utilities.testdata import get_local_testdata

import netCDF4 as nc4




CRS_LLDEG = 4326  # EPSG id of lat/lon (deg) coordinate referenence system (CRS)
CRS_CAEA = 3573  # EPSG id of equal-area    coordinate referenence system (CRS)
ROUTING_ID_FIELD = "HRU_ID"
NETCDF_INPUT_FIELD = "NetCDF_col"
DIM_NAMES = ("lon", "lat")
VAR_NAMES = ("lon", "lat")
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





# switched to local paths because with "get_local_testdata" following error:
#     RuntimeError: RAVENPY_TESTDATA_PATH env variable is not set
input_file = "/Users/j6mai/Documents/GitHub/GridWeightsGenerator/example/input_VIC/VIC_streaminputs.nc"
routing_file = "/Users/j6mai/Documents/GitHub/GridWeightsGenerator/example/maps/finalcat_hru_info.zip"

aggregate_forcings_to_HRUs( input_file, routing_file,
                                    dim_names=DIM_NAMES,
                                    var_names=VAR_NAMES,
                                    routing_id_field=ROUTING_ID_FIELD,
                                    netcdf_input_field=NETCDF_INPUT_FIELD,
                                    gauge_ids=None,
                                    sub_ids=None,
                                    area_error_threshold=AREA_ERROR_THRESHOLD,
                                    variables_to_aggregate=VARIABLES_TO_AGGREGATE )
