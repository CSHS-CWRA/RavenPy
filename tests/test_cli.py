from click.testing import CliRunner

from ravenpy.cli import generate_grid_weights

from .common import get_test_data


def test_generate_grid_weights():
    runner = CliRunner()
    ds = map(
        str,
        [
            get_test_data("raven-routing-sample", "finalcat_hru_info.zip")[0],
            get_test_data("raven-routing-sample", "VIC_streaminputs.nc")[0],
        ],
    )

    result = runner.invoke(generate_grid_weights, ds)

    assert result.exit_code == 0
    assert not result.exception
    assert ":NumberHRUs 51" in result.output
    assert ":NumberGridCells 100" in result.output
    assert len(result.output.split("\n")) == 216
