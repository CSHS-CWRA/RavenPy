import datetime as dt

import pytest
from pydantic import Field
from pydantic.dataclasses import dataclass

from ravenpy.new_config import commands as rc
from ravenpy.new_config import options as o
from ravenpy.new_config.base import Sym, SymConfig, Variable, encoder
from ravenpy.new_config.emulators.gr4jcn import GR4JCN
from ravenpy.new_config.rvs import Config

alt_names = {
    "PRECIP": "rain",
    "TEMP_MIN": "tmin",
    "TEMP_MAX": "tmax",
    "PET": "pet",
    "HYDROGRAPH": "qobs",
}


def test_emulator(dummy_config):
    cls, P = dummy_config
    conf = cls(Calendar="NOLEAP")

    assert conf.__config__.allow_mutation

    with pytest.raises(ValueError):
        assert conf.rvi

    # Set params
    num = conf.set_params([0.5])
    assert num.air_snow_coeff == 0.5

    with pytest.raises(ValueError):
        num.set_params([0.6])

    # Instantiate with numerical params
    nt = cls(params=[0.5], Calendar="NOLEAP")
    assert nt.air_snow_coeff == 0.5
