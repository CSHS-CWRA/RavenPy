from click.testing import CliRunner

from ravenpy.cli import generate_grid_weights

from .common import get_test_data


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
