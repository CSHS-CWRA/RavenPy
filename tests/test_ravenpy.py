from ravenpy import Emulator, parse, run


def test_run_and_parse(gr4jcn_config_rv):
    run(run_name="test", configdir=gr4jcn_config_rv, outputdir="output")
    out = parse(run_name="test", outputdir=gr4jcn_config_rv / "output")
    assert len(out["hydrograph"].q_sim.time) > 600


def test_emulator(gr4jcn_config, tmpdir):
    e = Emulator(config=gr4jcn_config, path=tmpdir)
    e.build()
    e.run()
    e.parse()

    assert ":RunName" in e.config.rvi
