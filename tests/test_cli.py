import re
from pathlib import Path

from click.testing import CliRunner

from ravenpy.cli import generate_grid_weights
from ravenpy.utilities.testdata import get_local_testdata


def test_generate_grid_weights_with_nc_input_and_2d_coords(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "bla.rvt"
    params = [
        get_local_testdata("raven-routing-sample/VIC_streaminputs.nc"),
        get_local_testdata("raven-routing-sample/finalcat_hru_info.zip"),
        "-o",
        output_path,
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception

    output = output_path.read_text()

    assert ":NumberHRUs 51" in output
    assert ":NumberGridCells 100" in output
    assert len(output.split("\n")) == 216

    # check derived weights
    weight = float(re.search("1 52 (.+)", output).group(1))
    assert abs(weight - 0.2610203097218425) < 1e-04


def test_generate_grid_weights_with_multiple_subids(tmp_path):
    # currently exactly same output as "test_generate_grid_weights_with_nc_input_and_2d_coords"
    # needs a "routing-file-path" with multiple gauges
    runner = CliRunner()
    output_path = tmp_path / "bla.rvt"
    params = [
        get_local_testdata("raven-routing-sample/VIC_streaminputs.nc"),
        get_local_testdata("raven-routing-sample/finalcat_hru_info.zip"),
        "-s",
        "7202",
        "-s",
        "6248",
        "-o",
        output_path,
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception

    output = output_path.read_text()

    assert ":NumberHRUs 51" in output
    assert ":NumberGridCells 100" in output
    assert len(output.split("\n")) == 216

    # check derived weights
    weight = float(re.search("1 52 (.+)", output).group(1))
    assert abs(weight - 0.2610203097218425) < 1e-04


def test_generate_grid_weights_with_nc_input_and_1d_coords(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "bla.rvt"
    params = [
        get_local_testdata("raven-routing-sample/era5-test-dataset-crop.nc"),
        get_local_testdata("raven-routing-sample/finalcat_hru_info.zip"),
        "--var-names",
        "longitude",
        "latitude",
        "-o",
        output_path,
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception

    output = output_path.read_text()

    assert ":NumberHRUs 51" in output
    assert ":NumberGridCells 9801" in output
    assert len(output.split("\n")) == 128

    # check derived weights
    weight = float(re.search("4 3731 (.+)", output).group(1))
    assert abs(weight - 0.0034512752779023515) < 1e-04


def test_generate_grid_weights_with_shp_input(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "bla.rvt"
    params = [
        get_local_testdata("raven-routing-sample/OTT_sub.zip"),
        get_local_testdata("raven-routing-sample/finalcat_hru_info.zip"),
        "-o",
        output_path,
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception

    output = output_path.read_text()

    assert ":NumberHRUs 51" in output
    assert ":NumberGridCells 810" in output
    assert len(output.split("\n")) == 230

    # check derived weights
    weight = float(re.search("13 238 (.+)", output).group(1))
    assert abs(weight - 0.5761414847779369) < 1e-04


def test_generate_grid_weights_with_weight_rescaling(tmp_path):
    runner = CliRunner()
    output_path = tmp_path / "bla.rvt"
    params = [
        get_local_testdata("raven-routing-sample/OTT_sub.zip"),
        get_local_testdata("raven-routing-sample/finalcat_hru_info.zip"),
        "--area-error-threshold",
        "0.42",
        "-o",
        output_path,
    ]
    params = map(str, params)

    result = runner.invoke(generate_grid_weights, params)

    assert result.exit_code == 0
    assert not result.exception

    output = output_path.read_text()

    assert ":NumberHRUs 51" in output
    assert ":NumberGridCells 810" in output
    assert len(output.split("\n")) == 230

    # check derived weights
    weight = float(re.search("13 238 (.+)", output).group(1))
    assert abs(weight - 0.9851111335377887) < 1e-04
