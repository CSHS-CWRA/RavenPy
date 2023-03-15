from typing import Literal, Tuple, Union

from ravenpy.config import options as o

from .commands import Process


# --- Processes --- #
class Precipitation(Process):
    """Precipitation"""

    algo: Literal["PRECIP_RAVEN", "RAVEN_DEFAULT"] = "PRECIP_RAVEN"


class CanopyEvap(Process):
    algo: Literal["CANEVP_RUTTER", "CANEVP_MAXIMUM"]


class SoilEvaporation(Process):
    """"""

    algo: Literal[
        "SOIL_EVAP_VIC",
        "SOIL_EVAP_HBV",
        "SOIL_EVAL_CHU",
        "SOILEVAP_GR4J",
        "SOIL_EVAP_TOPMODEL",
        "SOIL_EVAP_SEQUEN",
        "SOIL_EVAP_ROOTFRAC",
        "SOIL_EVAP_GOWSWER",
        "SOILEVAP_ALL",
    ]


class LakeEvaporation(Process):
    algo: Literal["LAKE_EVAP_BASIC"]


class OpenWaterEvaporation(Process):
    algo: Literal["OPEN_WATER_EVAP"]


class Infiltration(Process):
    """"""

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
    """"""

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
    """"""

    algo: Literal[
        "BASE_LINEAR",
        "BASE_POWER_LAW",
        "BASE_CONSTANT",
        "BASE_VIC",
        "BASE_THRESH_POWER",
        "BASE_GR4J",
        "BASE_TOPMODEL",
    ]


class Interflow(Process):
    algo: Literal["PRMS"]


class Seepage(Process):
    algo: Literal["SEEP_LINEAR"]


class DepressionOverflow(Process):
    algo: Literal["DFLOW_THRESHPOW", "DFLOW_LINEAR"]


class LakeRelease(Process):
    algo: Literal["LAKEREL_LINEAR"]


class Abstraction(Process):
    algo: Literal["ABST_PERCENTAGE", "ABST_FILL", "ABST_SCS"]


class SnowMelt(Process):
    algo: Literal["MELT_POTMELT"]


class SnowRefreeze(Process):
    algo: Literal["FREEZE_DEGREE_DAY"]


class SnowBalance(Process):
    algo: Literal[
        "SNOWBAL_SIMPLE_MELT",
        "SNOWBAL_COLD_CONTENT",
        "SNOWBAL_HBV",
        "SNOWBAL_TWO_LAYER",
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
    """"""

    algo: Literal["RAVEN_DEFAULT"]
    p: float = None


class Overflow(Process):
    algo: Literal["OVERFLOW_RAVEN"]
    _sub = "-->"
    _indent = "    "


class Split(Process):
    """"""

    p: float = None
    to: Tuple[str, str]


class Convolve(Process):
    """"""

    algo: Literal["CONVOL_GR4J_1", "CONVOL_GR4J_2", "CONVOL_GAMMA", "CONVOL_GAMMA_2"]


class LateralFlush(Process):
    """Lateral flush"""


class LateralEquilibrate(Process):
    """Lateral equilibrate

    Instantaneously equilibrates groundwater storage in basin HRUs.
    """


# --- End processes --- #
