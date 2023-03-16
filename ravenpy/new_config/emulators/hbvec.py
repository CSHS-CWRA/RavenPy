from typing import Dict, Literal, Sequence

from pydantic import Field
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable

import ravenpy.new_config.processes as p
from ravenpy.new_config import commands as rc
from ravenpy.new_config import options as o
from ravenpy.new_config.base import Params, Sym, SymConfig
from ravenpy.new_config.commands import HRU, PL, Process
from ravenpy.new_config.rvs import Config


@dataclass(config=SymConfig)
class P(Params):
    X01: Sym = Variable("X01")  #
    X02: Sym = Variable("X02")
    X03: Sym = Variable("X03")
    X04: Sym = Variable("X04")
    X05: Sym = Variable("X05")
    X06: Sym = Variable("X06")
    X07: Sym = Variable("X07")
    X08: Sym = Variable("X08")
    X09: Sym = Variable("X09")
    X10: Sym = Variable("X10")
    X11: Sym = Variable("X11")
    X12: Sym = Variable("X12")
    X13: Sym = Variable("X13")
    X14: Sym = Variable("X14")
    X15: Sym = Variable("X15")
    X16: Sym = Variable("X16")
    X17: Sym = Variable("X17")
    X18: Sym = Variable("X18")
    X19: Sym = Variable("X19")
    X20: Sym = Variable("X20")
    X21: Sym = Variable("X21")


class LandHRU(HRU):
    land_use_class: str = "LU_ALL"
    veg_class: str = "VEG_ALL"
    soil_profile: str = "DEFAULT_P"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    hru_type: Literal["land"] = "land"


class HRUs(rc.Command):
    """HRUs command for GR4J.

    Pydantic is able to automatically detect if an HRU is Land or Lake if `hru_type` is provided.
    """

    __root__: Sequence[LandHRU]


class HBVEC(Config):
    """Hydrologiska Byråns Vattenbalansavdelning – Environment Canada (HBV-EC)

    References
    ----------
    Lindström, G., et al., 1997. Development and test of the distributed HBV-96 hydrological model.
    Journal of Hydrology, 201 (1–4), 272–288. doi:10.1016/S0022-1694(97)00041-3

    Hamilton, A.S., Hutchinson, D.G., and Moore, R.D., 2000. Estimating winter streamflow using
    conceptual streamflow model. Journal of Cold Regions Engineering, 14 (4), 158–175.
    doi:10.1061/(ASCE)0887-381X(2000)14:4(158)

    Canadian Hydraulics Centre, 2010. Green kenue reference manual.
    Ottawa, Ontario: National Research Council.
    """

    params: P
    hrus: HRUs = Field(..., alias="HRUs")
    rain_snow_transition: rc.RainSnowTransition = Field(
        {"temp": P.X01, "delta": 2}, alias="RainSnowTransition"
    )
    global_parameter: Dict = Field(
        {
            "AdiabaticLapseRate": P.X13,
            "IrreducibleSnowSaturation": P.X04,
            "PRECIP_LAPSE": P.X12,
        },
        alias="GlobalParameters",
    )
    soil_classes: rc.SoilClasses = Field([])
