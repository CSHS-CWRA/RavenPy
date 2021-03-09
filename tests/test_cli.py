import re

import netCDF4 as nc4
from click.testing import CliRunner

from ravenpy.cli import aggregate_forcings_to_hrus, generate_grid_weights
from ravenpy.models.commands import GridWeightsCommand
from ravenpy.utilities.testdata import get_local_testdata


class TestGenerateGridWeights:
    def test_generate_grid_weights_with_nc_input_and_2d_coords(self, tmp_path):
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

    def test_generate_grid_weights_with_multiple_subids(self, tmp_path):
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

    def test_generate_grid_weights_with_nc_input_and_1d_coords(self, tmp_path):
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

    def test_generate_grid_weights_with_shp_input(self, tmp_path):
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

    def test_generate_grid_weights_with_weight_rescaling(self, tmp_path):
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


class TestAggregateForcingsToHRUs:
    def test_aggregate_forcings_to_hrus(self, tmp_path):
        runner = CliRunner()
        output_nc_file_path = tmp_path / "aggreg.nc"
        output_weight_file_path = tmp_path / "weight_aggreg.rvt"
        params = [
            get_local_testdata("raven-routing-sample/VIC_streaminputs.nc"),
            get_local_testdata("raven-routing-sample/VIC_streaminputs_weights.rvt"),
            "-v",
            "Streaminputs",
            "--output-nc-file",
            output_nc_file_path,
            "--output-weight-file",
            output_weight_file_path,
        ]
        params = map(str, params)

        result = runner.invoke(aggregate_forcings_to_hrus, params)

        assert result.exit_code == 0
        assert not result.exception

        output_rvt = output_weight_file_path.read_text()

        gws = GridWeightsCommand.parse(output_rvt)

        new_weights = gws.data

        # check new weights
        assert new_weights[0][0] == 1  # These are the HRU-IDs

        assert new_weights[0][1] == 0  # These need to be exactly [0,1,2,3,...,nHRU]
        assert new_weights[1][1] == 1  # These need to be exactly [0,1,2,3,...,nHRU]
        assert new_weights[2][1] == 2  # These need to be exactly [0,1,2,3,...,nHRU]
        assert new_weights[3][1] == 3  # These need to be exactly [0,1,2,3,...,nHRU]

        assert new_weights[0][2] == 1.0  # All new_weights[:][2] need to be 1.0
        assert new_weights[1][2] == 1.0  # All new_weights[:][2] need to be 1.0
        assert new_weights[2][2] == 1.0  # All new_weights[:][2] need to be 1.0
        assert new_weights[3][2] == 1.0  # All new_weights[:][2] need to be 1.0

        # check aggregated NetCDF file
        nc_in = nc4.Dataset(output_nc_file_path, "r")
        val = nc_in.variables["Streaminputs"][:]
        nc_in.close()

        assert abs(val[0, 0] - 0.017309) < 1e-04
        assert abs(val[16071, 0] - 0.569977) < 1e-04
        assert abs(val[0, 50] - 0.010276) < 1e-04
        assert abs(val[16071, 50] - 0.516639) < 1e-04
