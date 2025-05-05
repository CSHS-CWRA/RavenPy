from collections.abc import Sequence
from textwrap import dedent, indent
from typing import Literal, Optional

from .base import Sym
from .commands import Command, Process


# --- Processes --- #
class Precipitation(Process):
    """Precipitation."""

    algo: Literal["PRECIP_RAVEN", "RAVEN_DEFAULT"] = "PRECIP_RAVEN"


class CapillaryRise(Process):
    algo: Literal["RISE_HBV", "CRISE_HBV"]


class CanopyEvaporation(Process):
    algo: Literal["CANEVP_RUTTER", "CANEVP_MAXIMUM", "CANEVP_ALL"]


class CanopySublimation(Process):
    algo: Literal["CANEVP_ALL", "CANEVP_MAXIMUM", "CANSUBLIM_ALL", "CANSUBLIM_MAXIMUM"]


class SoilBalance(Process):
    algo: Literal["SOILBAL_SACSMA"]


class SoilEvaporation(Process):
    algo: Literal[
        "SOILEVAP_VIC",
        "SOILEVAP_HBV",
        "SOILEVAP_HYPR",
        "SOILEVAL_CHU",
        "SOILEVAP_UBC",
        "SOILEVAP_GR4J",
        "SOILEVAP_TOPMODEL",
        "SOILEVAP_SEQUEN",
        "SOILEVAP_ROOT",
        "SOILEVAP_ROOT_CONSTRAIN",
        "SOILEVAP_ROOTFRAC",
        "SOILEVAP_GAWSER",
        "SOILEVAP_FEDERER",
        "SOILEVAP_ALL",
        "SOILEVAP_LINEAR",
        "SOILEVAP_SACSMA",
        "SOILEVAP_HYMOD2",
    ]


class LakeEvaporation(Process):
    algo: Literal["LAKE_EVAP_BASIC"]


class LakeFreeze(Process):
    algo: Literal["LFREEZE_BASIC", "LFREEZE_THERMAL"]


class OpenWaterEvaporation(Process):
    algo: Literal["OPEN_WATER_EVAP", "OPEN_WATER_RIPARIAN"]


class Infiltration(Process):
    algo: Literal[
        "INF_RATIONAL",
        "INF_SCS",
        "INF_ALL_INFILTRATES",
        "INF_GR4J",
        "INF_GREEN_AMPT",
        "INF_GA_SIMPLE",
        "INF_UPSCALED_GREEN_AMPT",
        "INF_HBV",
        "INF_UBC",
        "INF_VIC",
        "INF_VIC_ARNO",
        "INF_PRMS",
        "INF_HMETS",
    ]


class Percolation(Process):
    algo: Literal[
        "PERC_GAWSER",
        "PERC_LINEAR",
        "PERC_POWER_LAW",
        "PERC_PRMS",
        "PERC_SACRAMENTO",
        "PERC_CONSTANT",
        "PERC_GR4J",
        "PERC_GR4JEXCH",
        "PERC_GR4JEXCH2",
    ]


class CapilaryRise(Process):
    algo: Literal["CRISE_HBV"]


class Baseflow(Process):
    algo: Literal[
        "BASE_LINEAR",
        "BASE_LINEAR_CONSTRAIN",
        "BASE_LINEAR_ANALYTIC",
        "BASE_POWER_LAW",
        "BASE_CONSTANT",
        "BASE_VIC",
        "BASE_THRESH_POWER",
        "BASE_THRESH_STOR",
        "BASE_GR4J",
        "BASE_TOPMODEL",
        "BASE_SACRAMENTO",
    ]


class Interflow(Process):
    algo: Literal["INTERFLOW_PRMS"]


class Recharge(Process):
    algo: Literal[
        "RECHARGE_FROMFILE",
        "RECHARGE_CONSTANT",
        "RECHARGE_DEFAULT",
        "RECHARGE_CONSTANT_OVERLAP",
        "RECHARGE_DATA",
        "RECHARGE_FLUX",
        "RAVEN_DEFAULT",
    ]


class Seepage(Process):
    algo: Literal["SEEP_LINEAR"]


class DepressionOverflow(Process):
    algo: Literal["DFLOW_THRESHPOW", "DFLOW_LINEAR"]


class LakeRelease(Process):
    algo: Literal["LAKEREL_LINEAR"]


class Abstraction(Process):
    algo: Literal["ABST_PERCENTAGE", "ABST_FILL", "ABST_SCS", "ABST_PDMROF"]


class SnowMelt(Process):
    algo: Literal["MELT_POTMELT"]


class SnowRefreeze(Process):
    algo: Literal["FREEZE_DEGREE_DAY"]


class SnowBalance(Process):
    algo: Literal[
        "SNOBAL_SIMPLE_MELT",
        "SNOBAL_COLD_CONTENT",
        "SNOBAL_HBV",
        "SNOBAL_TWO_LAYER",
        # "SNOWBAL_CEMA_NEIGE",
        "SNOBAL_CEMA_NIEGE",
        "SNOBAL_HMETS",
        "SNOWBAL_GAWSER",
        "SNOWBAL_UBC",
    ]


class Sublimation(Process):
    algo: Literal[
        "SUBLIM_SVERDRUP",
        "SUBLIM_KUZMIN",
        "SUBLIM_CENTRAL_SIERRA",
        "SUBLIM_PSBM",
        "SUBLIM_WILLIAMS",
    ]


class SnowAlbedoEvolve(Process):
    algo: Literal["SNOALB_UBC"]


class SnowTempEvolve(Process):
    algo: Literal["SNOTEMP_NEWTONS"]


class CanopyDrip(Process):
    algo: Literal["CANDRIP_RUTTER", "CANDRIP_SLOWDRAIN"]


class CropHeatUnitEvolve(Process):
    algo: Literal["CHU_ONTARIO"]


class GlacierMelt(Process):
    algo: Literal["GMELT_SIMPLE_MELT", "GMELT_HBV", "GMELT_UBC"]


class GlacierRelease(Process):
    algo: Literal["GRELEASE_LINEAR", "GRELEASE_HBV_EC"]


class Flush(Process):
    algo: Literal["RAVEN_DEFAULT"] = "RAVEN_DEFAULT"
    p: Optional[float] = None


class Overflow(Process):
    algo: Literal["OVERFLOW_RAVEN", "RAVEN_DEFAULT"]
    _sub = "-->"
    _indent = "    "


class Split(Process):
    p: float = None
    to: tuple[str, str]


class Convolve(Process):
    algo: Literal["CONVOL_GR4J_1", "CONVOL_GR4J_2", "CONVOL_GAMMA", "CONVOL_GAMMA_2"]


class LateralFlush(Process):
    """Lateral flush."""


class LateralEquilibrate(Process):
    """Lateral equilibrate.

    Instantaneously equilibrates groundwater storage in basin HRUs.
    """


# --- End processes --- #


class Conditional(Command):
    """Conditional statement."""

    kind: Literal["HRU_TYPE", "LAND_CLASS", "HRU_GROUP"]
    op: Literal["IS", "IS_NOT"]
    value: str

    _sub = "-->"

    def to_rv(self):
        cmd = self._sub + self.__class__.__name__
        return f":{cmd:<20} {self.kind} {self.op} {self.value}"


class ProcessGroup(Command):
    p: Sequence[Process]
    params: Sequence[Sym]

    @property
    def _template(self):
        return """
           :{_cmd}
           {_processes}
           :End{_cmd} CALCULATE_WTS {_params}
           """

    @property
    def _indent(self):
        return "    "

    def to_rv(self):
        d = {
            "_cmd": "ProcessGroup",
            "_processes": indent("\n".join(str(p) for p in self.p), self._indent),
            "_params": " ".join(map(str, self.params)),
        }
        return dedent(self._template).format(**d)
