import numpy as np
import pytest

from ravenpy.models.aggregate_forcings_to_HRUs import aggregate_forcings_to_HRUs
from ravenpy.utilities.testdata import get_local_testdata


def test_aggregate_forcings_to_HRUs():

    input_file = get_local_testdata("raven-routing-sample/VIC_streaminputs.nc")
    routing_file = get_local_testdata("raven-routing-sample/finalcat_hru_info.zip")

    DIM_NAMES = ("lon", "lat")
    VAR_NAMES = ("lon", "lat")
    ROUTING_ID_FIELD = "HRU_ID"
    NETCDF_INPUT_FIELD = "NetCDF_col"
    AREA_ERROR_THRESHOLD = 0.05
    VARIABLES_TO_AGGREGATE = ["Streaminputs"]

    agg_out = aggregate_forcings_to_HRUs( input_file, routing_file,
                                    dim_names=DIM_NAMES,
                                    var_names=VAR_NAMES,
                                    routing_id_field=ROUTING_ID_FIELD,
                                    netcdf_input_field=NETCDF_INPUT_FIELD,
                                    gauge_ids=None,
                                    sub_ids=None,
                                    area_error_threshold=AREA_ERROR_THRESHOLD,
                                    variables_to_aggregate=VARIABLES_TO_AGGREGATE )


    assert agg_out == 51   # dummy result for now :)

