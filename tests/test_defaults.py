from ravenpy.config import defaults


def test_defaults():
    import pint

    ureg = pint.UnitRegistry()

    for name, u in defaults.units.items():
        ureg(u)
