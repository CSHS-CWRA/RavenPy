import datetime as dt

import ravenpy.new_config.commands as rc
from ravenpy.new_config.rvs import RVC, RVH, RVI, RVP


class P:
    X01: float = 0
    X02: float = 0
    X03: float = 0
    X04: float = 0
    X05: float = 0
    X06: float = 0
    X07: float = 0
    X08: float = 0
    X09: float = 0
    X10: float = 0


rvp = RVP(
    SoilClasses=("TOPSOIL", "GWSOIL"),
    LandUseClasses=(
        rc.LandUseClass(name="LU_ALL", impermeable_frac=0, forest_coverage=1),
    ),
    VegetationClasses=(rc.VegetationClasses(name="VEG_ALL"),),
    SoilProfiles=(
        rc.SoilProfile(name="LAKE"),
        rc.SoilProfile(name="ROCK"),
        rc.SoilProfile(
            name="DEFAULT_P",
            soil_classes=("TOPSOIL", "GWSOIL"),
            thicknesses=(P.X05, 10.0),
        ),
    ),
    GlobalParameter={"TOC_MULTIPLIER": 1.0, "MOHYSE_PET_COEFFICIENT": P.X01},
    SoilParameterList=rc.SoilParameterList(
        names=[
            "POROSITY",
            "PET_CORRECTION",
            "HBV_BETA",
            "BASEFLOW_COEFF",
            "PERC_COEFF",
        ],
        records=[
            rc.PL(name="TOPSOIL", vals=[1, 1, 1, P.X07, P.X06]),
            rc.PL(name="GWSOIL", vals=[1, 1, 1, P.X08, 0]),
        ],
    ),
    LandUseParameterList=rc.LandUseParameterList(
        names=["MELT_FACTOR", "AET_COEFF", "FOREST_SPARSENESS", "DD_MELT_TEMP"],
        records=[rc.PL(name="DEFAULT", vals=[P.X03, P.X02, 0, P.X04])],
    ),
    VegetationParameterList=rc.VegetationParameterList(
        names=["SAI_HT_RATIO", "RAIN_ICEPT_PCT", "SNOW_ICEPT_PCT"],
        records=[rc.PL(name="[" "DEFAULT]", vals=[0, 0, 0])],
    ),
)


rvi = RVI(
    StartDate=dt.datetime(2000, 1, 1),
    SoilModel=2,
    PotentialMeltMethod="POTMELT_DEGREE_DAY",
    Routing="ROUTE_NONE",
    CatchmentRoute="ROUTE_GAMMA_CONVOLUTION",
    Evaporation="PET_MOHYSE",
    DirectEvaporation=True,
    RainSnowFraction="RAINSNOW_DATA",
)
