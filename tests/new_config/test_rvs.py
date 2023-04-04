import datetime as dt

import cftime
import pytest
from pydantic import ValidationError

from ravenpy.new_config.rvs import RVI, Config


def test_rvi_datetime():
    exp = cftime.datetime(1990, 1, 1, calendar="PROLEPTIC_GREGORIAN")

    rvi = RVI(start_date=dt.datetime(1990, 1, 1))
    assert rvi.start_date == exp

    rvi = RVI(start_date="1990-01-01")
    assert rvi.start_date == exp

    rvi = RVI(start_date="1990-01-01T00:00:00")
    assert rvi.start_date == exp

    rvi = RVI(end_date=cftime.datetime(1990, 1, 1))
    assert rvi.end_date == exp

    with pytest.raises(ValidationError):
        rvi = RVI(start_date=(dt.datetime(1990, 1, 1),))


def test_duplicate():
    conf = Config(start_date="1990-01-01")

    # Updating values with an alias and an attribute name
    out = conf.duplicate(Duration=10, debug_mode=True)

    assert isinstance(out.start_date, cftime.datetime)
    assert out.duration == 10
    assert out.debug_mode
