from ravenpy.config import defaults


def test_defaults():
    from pint import unit

    for name, u in defaults.units.items():
        unit.Unit(u)
