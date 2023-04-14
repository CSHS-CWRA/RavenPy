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

    # Create list of output paths using glob
    paths = p.glob("**/output")

    ens = EnsembleReader(paths=paths, dim="parameters")
    assert len(ens.hydrograph.parameters) == 2
