import pathlib
from importlib.util import find_spec

from ravenpy import Emulator, EnsembleReader


def test_ensemble_reader(gr4jcn_config, tmp_path):
    """Mimics code in docs/outputs.md"""
    conf, params = gr4jcn_config

    params = [params, params]

    # Output directory for all simulations
    p = tmp_path / "ensemble"

    # Writing to the same custom output cause segfaults
    # conf.custom_output = [rc.CustomOutput(time_per="DAILY", stat="AVERAGE", variable="RAINFALL",
    # space_agg="ENTIRE_WATERSHED")]

    # Run the model for each parameter set in `params`
    runs = [
        Emulator(conf.set_params(param), workdir=p / f"m{i}").run()
        for i, param in enumerate(params)
    ]

    ens = EnsembleReader(runs=runs, dim="parameters")
    assert len(ens.hydrograph.parameters) == 2

    # Create a list of output paths using glob
    paths = p.glob("**/output")

    ens = EnsembleReader(paths=paths, dim="parameters")
    assert len(ens.hydrograph.parameters) == 2


def test_package_metadata():
    """Test the package metadata."""
    project = find_spec("ravenpy")

    assert project is not None
    assert project.submodule_search_locations is not None

    location = project.submodule_search_locations[0]

    metadata = pathlib.Path(location).resolve().joinpath("__init__.py")

    with metadata.open() as f:
        contents = f.read()
        assert """David Huard""" in contents
        assert '__email__ = "huard.david@ouranos.ca"' in contents
        assert '__version__ = "0.19.0"' in contents
