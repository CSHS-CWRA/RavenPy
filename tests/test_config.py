from ravenpy.config import defaults


def test_gr4jcn(tmpdir):
    from ravenpy.models.emulators import GR4JCN
    m = GR4JCN(params=[1,2,3,4,5,6])
    m.write(tmpdir)


def test_defaults():
    import pint

    ureg = pint.UnitRegistry()

    for name, u in defaults.units.items():
        ureg(u)
