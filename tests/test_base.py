from ravenpy.config import commands as rc


def test_line_command():
    c = rc.SBGroupPropertyMultiplier(
        group_name="Land", parameter_name="MANNINGS_N", mult=1.0
    )

    out = c.to_rv()
    assert out == ":SBGroupPropertyMultiplier Land MANNINGS_N 1.0\n"
