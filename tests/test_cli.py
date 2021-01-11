import re

from click.testing import CliRunner

from ravenpy.cli import generate_grid_weights
from ravenpy.utilities.testdata import get_test_data


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
    weight = float(re.search("1 52 (.+)", result.output).group(1))
    assert abs(weight - 0.2610203097218425) < 1e-04


def test_generate_grid_weights_with_multiple_subids():
    # currently exactly same output as "test_generate_grid_weights_with_nc_input_and_2d_coords"
    # needs a "routing-file-path" with multiple gauges
    runner = CliRunner()
    params = [
        get_test_data("raven-routing-sample", "VIC_streaminputs.nc")[0],
        get_test_data("raven-routing-sample", "finalcat_hru_info.zip")[0],
        "-s",
        "7202",
        "-s",
        "6248",
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception
    assert ":NumberHRUs 51" in result.output
    assert ":NumberGridCells 100" in result.output
    assert len(result.output.split("\n")) == 216

    # check derived weights
    weight = float(re.search("1 52 (.+)", result.output).group(1))
    assert abs(weight - 0.2610203097218425) < 1e-04


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
    weight = float(re.search("4 3731 (.+)", result.output).group(1))
    assert abs(weight - 0.0034512752779023515) < 1e-04


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
    weight = float(re.search("13 238 (.+)", result.output).group(1))
    assert abs(weight - 0.5761414847779369) < 1e-04


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
    weight = float(re.search("13 238 (.+)", result.output).group(1))
    assert abs(weight - 0.9851111335377887) < 1e-04
