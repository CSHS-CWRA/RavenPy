from ravenpy.config import defaults


def test_defaults():
    import pint

    ureg = pint.UnitRegistry()

    for u in defaults.units.values():
        ureg(u)
