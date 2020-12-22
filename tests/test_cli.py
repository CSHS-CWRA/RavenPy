from click.testing import CliRunner

from ravenpy.cli import generate_grid_weights

from .common import get_test_data

import numpy as np


def test_generate_grid_weights_with_nc_input_and_2d_coords():
    runner = CliRunner()
    params = [
        get_test_data("raven-routing-sample", "VIC_streaminputs.nc")[0],
        get_test_data("raven-routing-sample", "finalcat_hru_info.zip")[0],
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception
    assert ":NumberHRUs 51" in result.output
    assert ":NumberGridCells 100" in result.output   
    assert len(result.output.split("\n")) == 216

    # check derived weights
    target_HRU_ID = 1
    target_cell_ID = 52
    target_weight = 0.2610203097218425
    result_as_array = result.output.split("\n")[3:-2]  # remove header and footer lines 
    result_as_array = np.array( [  list(map(float,ii.strip().split())) for ii in result_as_array ] ) # convert to array
    weight = result_as_array[np.logical_and(result_as_array[:,0]==target_HRU_ID, result_as_array[:,1]==target_cell_ID),2][0] # find row with predefined HRU_ID and Cell_ID and extract weight
    assert np.isclose(weight, target_weight, rtol=1e-04, atol=1e-04, equal_nan=False)


def test_generate_grid_weights_with_multiple_subids():
    # currently exactly same output as "test_generate_grid_weights_with_nc_input_and_2d_coords"
    # needs a "routing-file-path" with multiple gauges
    runner = CliRunner()
    params = [
        get_test_data("raven-routing-sample", "VIC_streaminputs.nc")[0],
        get_test_data("raven-routing-sample", "finalcat_hru_info.zip")[0],
        "--sub-id",
        "7202,6248",
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception
    assert ":NumberHRUs 51" in result.output
    assert ":NumberGridCells 100" in result.output    
    assert len(result.output.split("\n")) == 216

    # check derived weights
    target_HRU_ID = 1
    target_cell_ID = 52
    target_weight = 0.2610203097218425
    result_as_array = result.output.split("\n")[3:-2]  # remove header and footer lines 
    result_as_array = np.array( [  list(map(float,ii.strip().split())) for ii in result_as_array ] ) # convert to array
    weight = result_as_array[np.logical_and(result_as_array[:,0]==target_HRU_ID, result_as_array[:,1]==target_cell_ID),2][0] # find row with predefined HRU_ID and Cell_ID and extract weight
    assert np.isclose(weight, target_weight, rtol=1e-04, atol=1e-04, equal_nan=False)


def test_generate_grid_weights_with_nc_input_and_1d_coords():
    runner = CliRunner()
    params = [
        get_test_data("raven-routing-sample", "era5-test-dataset-crop.nc")[0],
        get_test_data("raven-routing-sample", "finalcat_hru_info.zip")[0],
        "--var-names",
        "longitude",
        "latitude",
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception
    assert ":NumberHRUs 51" in result.output
    assert ":NumberGridCells 9801" in result.output
    assert len(result.output.split("\n")) == 128

    # check derived weights
    target_HRU_ID = 4
    target_cell_ID = 3731
    target_weight = 0.0034512752779023515
    result_as_array = result.output.split("\n")[3:-2]  # remove header and footer lines 
    result_as_array = np.array( [  list(map(float,ii.strip().split())) for ii in result_as_array ] ) # convert to array
    weight = result_as_array[np.logical_and(result_as_array[:,0]==target_HRU_ID, result_as_array[:,1]==target_cell_ID),2][0] # find row with predefined HRU_ID and Cell_ID and extract weight
    assert np.isclose(weight, target_weight, rtol=1e-04, atol=1e-04, equal_nan=False)


def test_generate_grid_weights_with_shp_input():
    runner = CliRunner()
    params = [
        get_test_data("raven-routing-sample", "OTT_sub.zip")[0],
        get_test_data("raven-routing-sample", "finalcat_hru_info.zip")[0],
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception
    assert ":NumberHRUs 51" in result.output
    assert ":NumberGridCells 810" in result.output
    assert len(result.output.split("\n")) == 230

    # check derived weights
    target_HRU_ID = 13
    target_cell_ID = 238
    target_weight = 0.5761414847779369
    result_as_array = result.output.split("\n")[3:-2]  # remove header and footer lines 
    result_as_array = np.array( [  list(map(float,ii.strip().split())) for ii in result_as_array ] ) # convert to array
    weight = result_as_array[np.logical_and(result_as_array[:,0]==target_HRU_ID, result_as_array[:,1]==target_cell_ID),2][0] # find row with predefined HRU_ID and Cell_ID and extract weight
    assert np.isclose(weight, target_weight, rtol=1e-04, atol=1e-04, equal_nan=False)


def test_generate_grid_weights_with_weight_rescaling():
    runner = CliRunner()
    params = [
        get_test_data("raven-routing-sample", "OTT_sub.zip")[0],
        get_test_data("raven-routing-sample", "finalcat_hru_info.zip")[0],
        "--area-error-threshold",
        "0.42",
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception
    assert ":NumberHRUs 51" in result.output
    assert ":NumberGridCells 810" in result.output
    assert len(result.output.split("\n")) == 230

    # check derived weights
    target_HRU_ID = 13
    target_cell_ID = 238
    target_weight = 0.9851111335377887
    result_as_array = result.output.split("\n")[3:-2]  # remove header and footer lines 
    result_as_array = np.array( [  list(map(float,ii.strip().split())) for ii in result_as_array ] ) # convert to array
    weight = result_as_array[np.logical_and(result_as_array[:,0]==target_HRU_ID, result_as_array[:,1]==target_cell_ID),2][0] # find row with predefined HRU_ID and Cell_ID and extract weight
    assert np.isclose(weight, target_weight, rtol=1e-04, atol=1e-04, equal_nan=False)
