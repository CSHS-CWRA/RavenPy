from ravenpy import parse, run


def test_run_and_parse(gr4jcn_config):
    run(identifier="test", configdir=gr4jcn_config, outputdir="output")
    out = parse(identifier="test", outputdir=gr4jcn_config / "output")
    assert len(out["hydrograph"].q_sim.time) > 600
