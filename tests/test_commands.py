import pathlib
import re
from collections.abc import Sequence
from textwrap import dedent
from typing import Union

import pytest
import xarray as xr
from pydantic import Field

from ravenpy.config import commands as rc
from ravenpy.config import options as o
from ravenpy.config.base import RV, optfield

salmon_land_hru_1 = dict(
    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659, hru_type="land"
)


def test_custom_output():
    class Test(RV):
        co: Sequence[rc.CustomOutput] = optfield(alias="CustomOutput")

    co = rc.CustomOutput(
        time_per="YEARLY",
        stat="AVERAGE",
        variable="PRECIP",
        space_agg="ENTIRE_WATERSHED",
    )
    assert (
        Test(co=[co]).to_rv().strip()
        == ":CustomOutput         YEARLY AVERAGE PRECIP ENTIRE_WATERSHED"
    )

    new = rc.CustomOutput.parse(str(co))
    assert new == co


def test_evaluation_metrics():

    evaluation_metrics = [
        "ABSERR",
        "ABSERR_RUN",
        "ABSMAX",
        "DIAG_SPEARMAN",
        "FUZZY_NASH",
        "KGE_PRIME",
        "KLING_GUPTA",
        "LOG_NASH",
        "NASH_SUTCLIFFE",
        "NSC",
        "PCT_BIAS",
        "PDIFF",
        "RCOEFF",
        "RMSE",
        "TMVOL",
    ]

    class Test(RV):
        em: Sequence[o.EvaluationMetrics] = Field(
            evaluation_metrics, alias="EvaluationMetrics"
        )

    assert (
        Test().to_rv().strip()
        == f":EvaluationMetrics    {' '.join(evaluation_metrics)}"
    )


def test_hrus():
    hrus = rc.HRUs([salmon_land_hru_1])
    hrus.to_rv()
    assert hrus[0].subbasin_id == 1


def test_hrus_parse(yangtze):
    f = yangtze.fetch("clrh/mattawin/06FB002.rvh")
    hrus = rc.HRUs.parse(pathlib.Path(f).read_text())
    assert len(hrus) == 40
    hru = hrus[0]

    assert hru.hru_id == 38
    assert hru.area == 2.1924
    assert hru.elevation == 1440.5898
    assert hru.latitude == 46.6007
    assert hru.longitude == -73.9222
    assert hru.subbasin_id == 23007946
    assert hru.land_use_class == "LAKE"
    assert hru.veg_class == "LAKE"
    assert hru.soil_profile == "LAKE"
    assert hru.aquifer_profile == "[NONE]"
    assert hru.terrain_class == "[NONE]"
    assert hru.slope == 30.986380
    assert hru.aspect == 190.30589321792877


def test_hru_state():
    s = rc.HRUState(hru_id=1, data={"SOIL[0]": 1, "SOIL[1]": 2.0})
    assert str(s) == "1,1.0,2.0"


def test_linear_transform():
    assert rc.LinearTransform(scale=1, offset=0).to_rv() == ""


def test_global_parameter():
    class Test(RV):
        gp: dict = Field({"A": 1}, alias="GlobalParameter")

    s = Test().to_rv().strip()

    assert s == ":GlobalParameter A 1"


def test_vegetation_classes():
    class Test(RV):
        vegetation_classes: rc.VegetationClasses = Field(
            rc.VegetationClasses([{"name": "VEG_ALL"}, {"name": "VEG_WATER"}]),
            alias="VegetationClasses",
        )

    s = dedent(Test().to_rv().strip())
    assert (
        s
        == dedent(
            """
    :VegetationClasses
      :Attributes     ,            MAX_HT,           MAX_LAI,     MAX_LEAF_COND
      :Units          ,                 m,              none,          mm_per_s
      VEG_ALL         ,               0.0,               0.0,               0.0
      VEG_WATER       ,               0.0,               0.0,               0.0
    :EndVegetationClasses"""
        ).strip()
    )


def test_land_use_class():
    class Test(RV):
        lu_classes: rc.LandUseClasses = Field(None, alias="LandUseClasses")

    t = Test(
        LandUseClasses=[
            {"name": "Lake"},
            {"name": "Land", "impermeable_frac": 0.1, "forest_coverage": 0.8},
        ]
    )
    s = dedent(t.to_rv()).strip()
    assert (
        s
        == dedent(
            """
        :LandUseClasses
          :Attributes     ,  IMPERMEABLE_FRAC,   FOREST_COVERAGE
          :Units          ,             fract,             fract
          Lake            ,               0.0,               0.0
          Land            ,               0.1,               0.8
        :EndLandUseClasses
        """
        ).strip()
    )


def test_soil_classes():
    from ravenpy.config.commands import SoilClasses

    c = SoilClasses([{"name": "TOPSOIL"}, {"name": "FAST_RES"}, {"name": "SLOW_RES"}])
    assert dedent(c.to_rv()) == dedent(
        """
    :SoilClasses
      :Attributes     ,             %SAND,             %CLAY,             %SILT,          %ORGANIC
      :Units          ,              none,              none,              none,              none
      TOPSOIL
      FAST_RES
      SLOW_RES
    :EndSoilClasses
    """
    )

    c = SoilClasses(
        [
            {"name": "TOPSOIL", "mineral": (0.1, 0.7, 0.2), "organic": 0.1},
        ]
    )
    assert dedent(c.to_rv()) == dedent(
        """
    :SoilClasses
      :Attributes     ,             %SAND,             %CLAY,             %SILT,          %ORGANIC
      :Units          ,              none,              none,              none,              none
      TOPSOIL         ,               0.1,               0.7,               0.2,               0.1
    :EndSoilClasses
    """
    )


def test_soil_model():
    class Test(RV):
        soil_model: rc.SoilModel = Field(None, alias="SoilModel")

    t = Test(soil_model=3)
    assert t.to_rv().strip() == ":SoilModel            SOIL_MULTILAYER 3"


def test_soil_profiles(x=42.0):
    c = rc.SoilProfiles(
        [
            {
                "name": "DEFAULT_P",
                "soil_classes": ["TOPSOIL", "FAST_RES", "SLOW_RES"],
                "thicknesses": [x, 100.0, 100.0],
            },
        ]
    )

    pat = r"DEFAULT_P\s*,\s*3,\s*TOPSOIL,\s*(.+),\s*FAST_RES,\s*100.0,\s*SLOW_RES,\s*100.0"

    m = re.findall(pat, str(c))[0]
    assert m == "42.0"


def test_hydrologic_processes():
    from ravenpy.config.processes import Precipitation, SnowTempEvolve

    hp = [
        Precipitation(algo="PRECIP_RAVEN", source="ATMOS_PRECIP", to="MULTIPLE"),
        SnowTempEvolve(algo="SNOTEMP_NEWTONS", source="SNOW_TEMP"),
    ]

    class Test(RV):
        hydrologic_processes: Sequence[rc.Process] = Field(
            hp, alias="HydrologicProcesses"
        )

    out = dedent(Test().to_rv().strip())
    assert (
        out
        == dedent(
            """
        :HydrologicProcesses
          :Precipitation        PRECIP_RAVEN        ATMOS_PRECIP        MULTIPLE
          :SnowTempEvolve       SNOTEMP_NEWTONS     SNOW_TEMP
        :EndHydrologicProcesses
        """
        ).strip()
    )


def test_process_group():
    from ravenpy.config.processes import Precipitation, ProcessGroup, SnowTempEvolve

    hp = [
        Precipitation(algo="PRECIP_RAVEN", source="ATMOS_PRECIP", to="MULTIPLE"),
        SnowTempEvolve(algo="SNOTEMP_NEWTONS", source="SNOW_TEMP"),
    ]
    c = ProcessGroup(p=hp, params=(0, 1))
    out = dedent(c.to_rv().strip())

    assert (
        out
        == dedent(
            """
    :ProcessGroup
        :Precipitation        PRECIP_RAVEN        ATMOS_PRECIP        MULTIPLE
        :SnowTempEvolve       SNOTEMP_NEWTONS     SNOW_TEMP
    :EndProcessGroup CALCULATE_WTS 0.0 1.0
        """
        ).strip()
    )


def test_soil_parameter_list():
    c = rc.SoilParameterList(
        parameters=["POROSITY", "FIELD_CAPACITY"],
        pl=[
            rc.PL(name="[DEFAULT]", values=[1, 0]),
            rc.PL(name="FAST_RES", values=[None, None]),
        ],
    )

    assert dedent(c.to_rv()) == dedent(
        """
    :SoilParameterList
      :Parameters     ,          POROSITY,    FIELD_CAPACITY
      :Units          ,              none,              none
      [DEFAULT]       ,               1.0,               0.0
      FAST_RES        ,          _DEFAULT,          _DEFAULT
    :EndSoilParameterList
    """
    )


def test_vegetation_parameter_list():
    c = rc.VegetationParameterList(
        parameters=[
            "MAX_CAPACITY",
            "MAX_SNOW_CAPACITY",
            "RAIN_ICEPT_PCT",
            "SNOW_ICEPT_PCT",
        ],
        pl=[rc.PL(name="VEG_ALL", values=[10000, 10000, 0.88, 0.88])],
    )

    assert dedent(c.to_rv()) == dedent(
        """
    :VegetationParameterList
      :Parameters     ,      MAX_CAPACITY, MAX_SNOW_CAPACITY,    RAIN_ICEPT_PCT,    SNOW_ICEPT_PCT
      :Units          ,              none,              none,              none,              none
      VEG_ALL         ,           10000.0,           10000.0,              0.88,              0.88
    :EndVegetationParameterList
    """
    )


def test_land_use_parameter_list():
    c = rc.LandUseParameterList(
        parameters=["MELT_FACTOR", "MIN_MELT_FACTOR"],
        pl=[rc.PL(name="[DEFAULT]", values=[1, 2.2])],
    )
    assert dedent(c.to_rv()) == dedent(
        """
    :LandUseParameterList
      :Parameters     ,       MELT_FACTOR,   MIN_MELT_FACTOR
      :Units          ,              none,              none
      [DEFAULT]       ,               1.0,               2.2
    :EndLandUseParameterList
    """
    )

    with pytest.raises(ValueError):
        rc.LandUseParameterList(
            names=["MELT_FACTOR", "MIN_MELT_FACTOR"],
            pl=[rc.PL(name="[DEFAULT]", values=[None, 2.2])],
        )


class TestHRUStateVariableTable:
    def test_to_rv(self):
        s1 = rc.HRUState(hru_id=1, data={"SOIL[0]": 0.1, "SOIL[1]": 1.0})
        s2 = rc.HRUState(hru_id=2, data={"SOIL[3]": 3, "SOIL[2]": 2.0})
        t = rc.HRUStateVariableTable([s1, s2])
        assert dedent(t.to_rv()) == dedent(
            """
            :HRUStateVariableTable
              :Attributes     ,           SOIL[0],           SOIL[1],           SOIL[2],           SOIL[3]
              1,0.1,1.0,0.0,0.0
              2,0.0,0.0,2.0,3.0
            :EndHRUStateVariableTable
            """
        )

        s1 = dict(hru_id=1, data={"SOIL[0]": 0.1, "SOIL[1]": 1.0})
        rc.HRUStateVariableTable([s1])

    def test_parse(self):
        solution = """
:HRUStateVariableTable
  :Attributes,SURFACE_WATER,ATMOSPHERE,ATMOS_PRECIP,PONDED_WATER,SOIL[0],SOIL[1],SNOW,SNOW_LIQ,CUM_SNOWMELT,CONVOLUTION[0],CONVOLUTION[1],AET,CONV_STOR[0],CONV_STOR[1],CONV_STOR[2],CONV_STOR[3],CONV_STOR[4],CONV_STOR[5],CONV_STOR[6],CONV_STOR[7],CONV_STOR[8],CONV_STOR[9],CONV_STOR[10],CONV_STOR[11],CONV_STOR[12],CONV_STOR[13],CONV_STOR[14],CONV_STOR[15],CONV_STOR[16],CONV_STOR[17],CONV_STOR[18],CONV_STOR[19],CONV_STOR[20],CONV_STOR[21],CONV_STOR[22],CONV_STOR[23],CONV_STOR[24],CONV_STOR[25],CONV_STOR[26],CONV_STOR[27],CONV_STOR[28],CONV_STOR[29],CONV_STOR[30],CONV_STOR[31],CONV_STOR[32],CONV_STOR[33],CONV_STOR[34],CONV_STOR[35],CONV_STOR[36],CONV_STOR[37],CONV_STOR[38],CONV_STOR[39],CONV_STOR[40],CONV_STOR[41],CONV_STOR[42],CONV_STOR[43],CONV_STOR[44],CONV_STOR[45],CONV_STOR[46],CONV_STOR[47],CONV_STOR[48],CONV_STOR[49],CONV_STOR[50],CONV_STOR[51],CONV_STOR[52],CONV_STOR[53],CONV_STOR[54],CONV_STOR[55],CONV_STOR[56],CONV_STOR[57],CONV_STOR[58],CONV_STOR[59],CONV_STOR[60],CONV_STOR[61],CONV_STOR[62],CONV_STOR[63],CONV_STOR[64],CONV_STOR[65],CONV_STOR[66],CONV_STOR[67],CONV_STOR[68],CONV_STOR[69],CONV_STOR[70],CONV_STOR[71],CONV_STOR[72],CONV_STOR[73],CONV_STOR[74],CONV_STOR[75],CONV_STOR[76],CONV_STOR[77],CONV_STOR[78],CONV_STOR[79],CONV_STOR[80],CONV_STOR[81],CONV_STOR[82],CONV_STOR[83],CONV_STOR[84],CONV_STOR[85],CONV_STOR[86],CONV_STOR[87],CONV_STOR[88],CONV_STOR[89],CONV_STOR[90],CONV_STOR[91],CONV_STOR[92],CONV_STOR[93],CONV_STOR[94],CONV_STOR[95],CONV_STOR[96],CONV_STOR[97],CONV_STOR[98],CONV_STOR[99]
  :Units,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm
  1,0.00000,0.00000,-0.16005,0.00000,144.54883,455.15708,0.16005,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000
:EndHRUStateVariableTable
        """
        sv = rc.HRUStateVariableTable.parse(solution).root
        assert len(sv) == 1
        assert sv[0].hru_id == 1
        assert sv[0].data["ATMOS_PRECIP"] == -0.16005


def test_read_from_netcdf(yangtze):
    from ravenpy.config.commands import ReadFromNetCDF

    f = yangtze.fetch(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    c = ReadFromNetCDF.from_nc(f, "PRECIP", station_idx=1, alt_names=("rain",))
    s = dedent(c.to_rv())

    pat = re.compile(
        r"""
    :ReadFromNetCDF\s+
      :FileNameNC\s+.+/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc\s+
      :VarNameNC\s+rain\s+
      :DimNamesNC\s+time\s+
      :StationIdx\s+1\s+
      :LatitudeVarNameNC\s+lat\s+
      :LongitudeVarNameNC\s+lon\s+
    :EndReadFromNetCDF
    """,
        re.VERBOSE + re.MULTILINE,
    )
    assert pat.search(s) is not None
    assert isinstance(c.da, xr.DataArray)


def test_station_forcing(yangtze):
    f = yangtze.fetch(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc"
    )
    c = rc.StationForcing.from_nc(f, "PRECIP", station_idx=1, alt_names="rain")
    dedent(c.to_rv())


def test_gridded_forcing(yangtze):
    # TODO: Make sure dimensions are in the order x, y, t.
    fn = yangtze.fetch("raven-routing-sample/VIC_temperatures.nc")
    rc.GriddedForcing.from_nc(fn, data_type="TEMP_AVE", alt_names=("Avg_temp",))
    # assert gf.dim_names_nc == ("lon_dim", "lat_dim", "time")

    fn = yangtze.fetch("raven-routing-sample/VIC_streaminputs.nc")
    rc.GriddedForcing.from_nc(fn, data_type="PRECIP", alt_names=("Streaminputs",))
    # assert gf.dim_names_nc == ("lon_dim", "lat_dim", "time")

    fn = yangtze.fetch("cmip5/tas_Amon_CanESM2_rcp85_r1i1p1_200601-210012_subset.nc")
    gf = rc.GriddedForcing.from_nc(fn, data_type="TEMP_AVE", engine="netcdf4")
    assert gf.latitude_var_name_nc is None


def test_gauge(yangtze, tmp_path):
    salmon_file = yangtze.fetch(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    salmon_file_tmp = tmp_path / "salmon_river_near_prince_george-tmp.nc"
    salmon_file_tmp.write_bytes(pathlib.Path(salmon_file).read_bytes())

    g = rc.Gauge.from_nc(
        salmon_file_tmp,
        alt_names={"RAINFALL": "rain", "SNOWFALL": "snow"},
        data_kwds={"ALL": {"Deaccumulate": True}},
    )

    assert "Data" in g.to_rv()
    assert isinstance(g.ds, xr.Dataset)
    assert "RAINFALL" in g.ds
    assert g.data[0].read_from_netcdf.deaccumulate


def test_gauge_raises(yangtze):
    f = yangtze.fetch(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )

    with pytest.raises(ValueError):
        rc.Gauge.from_nc(
            f,
            data_type=[
                "RAINFALL",
            ],
            alt_names={"RAINFALL": "bad_name"},
        )


def test_grid_weights():
    gw = rc.GridWeights()
    txt = gw.to_rv()
    assert len(txt.splitlines()) == 6

    parsed = rc.GridWeights.parse(txt)
    assert parsed.number_hrus == gw.number_hrus
    assert parsed.number_grid_cells == gw.number_grid_cells
    assert parsed.data == gw.data


def test_redirect_to_file(yangtze):
    f = yangtze.fetch("raven-routing-sample/VIC_test_nodata_weights.rvt")
    r = rc.RedirectToFile(f)
    assert re.match(r"^:RedirectToFile\s+(\S+)$", r.to_rv())

    class Test(RV):
        gw: Union[rc.GridWeights, rc.RedirectToFile] = Field(
            None, alias="RedirectToFile"
        )

    t = Test(gw=f)
    assert t.gw.root == pathlib.Path(f)


def test_subbasin_properties():
    c = rc.SubBasinProperties(
        parameters=["GAMMA_SCALE", "GAMMA_SHAPE"],
        records=[{"sb_id": 1, "values": (0, 1)}],
    )
    assert "0, 1" in c.to_rv()


def test_subbasins():
    c = rc.SubBasins([{"name": "SB1"}])
    c.to_rv()
    assert c[0].gauged

    class Test(RV):
        sub_basins: rc.SubBasins = Field([rc.SubBasin()], alias="SubBasins")
        hrus: rc.HRUs = Field([rc.HRU()], alias="HRUs")

    Test().to_rv()


def test_subbasins_parse(yangtze):
    f = yangtze.fetch("clrh/mattawin/06FB002.rvh")
    sb = rc.SubBasins.parse(pathlib.Path(f).read_text())
    assert len(sb) == 35


def test_reservoir_parse(yangtze):
    f = yangtze.fetch("clrh/mattawin/Lakes.rvh")
    rs = rc.Reservoir.parse(pathlib.Path(f).read_text())
    assert len(rs) == 5
    r = rs[0]

    assert r.name == "Lake_107056"
    assert r.subbasin_id == 23007946
    assert r.hru_id == 38
    assert r.type == "RESROUTE_STANDARD"
    assert r.weir_coefficient == 0.6
    assert r.crest_width == 7.7399
    assert r.max_depth == 9.6
    assert r.lake_area == 2192395.4127609455
    assert r.seepage_parameters.k_seep == 0
    assert r.seepage_parameters.h_ref == 0


def test_ensemble_mode():
    c = rc.EnsembleMode(n=10)
    s = c.to_rv().strip()
    assert s == ":EnsembleMode         ENSEMBLE_ENKF 10"


def test_observation_error_model():
    c = rc.ObservationErrorModel(
        state="STREAMFLOW", dist="DIST_UNIFORM", p1=1, p2=2, adj="ADDITIVE"
    )

    assert (
        c.to_rv().strip()
        == ":ObservationErrorModel STREAMFLOW DIST_UNIFORM 1.0 2.0 ADDITIVE"
    )


def test_forcing_perturbation():
    c = rc.ForcingPerturbation(
        forcing="TEMP_AVE", dist="DIST_UNIFORM", p1=1, p2=2, adj="ADDITIVE"
    )

    assert (
        c.to_rv().strip()
        == ":ForcingPerturbation  TEMP_AVE DIST_UNIFORM 1.0 2.0 ADDITIVE"
    )


def test_hru_group():
    c = rc.HRUGroup(name="ForestedHRUs", groups=["1", "3-5"])
    s = dedent(c.to_rv()).strip()
    assert (
        s
        == dedent(
            """
            :HRUGroup ForestedHRUs
              1,3-5
            :EndHRUGroup
                """
        ).strip()
    )


def test_subbasin_group_parse(yangtze):
    f = yangtze.fetch("clrh/mattawin/06FB002.rvh")
    sbgs = rc.SubBasinGroup.parse(pathlib.Path(f).read_text())
    assert len(sbgs) == 2
    sbg = sbgs[0]
    assert sbg.name == "Allsubbasins"
    assert len(sbg.sb_ids) == 35
    sbg.sb_ids[0] == 23007946


def test_channel_profile_parse(yangtze):
    f = yangtze.fetch("clrh/mattawin/channel_properties.rvp")
    cps = rc.ChannelProfile.parse(pathlib.Path(f).read_text())
    assert len(cps) == 20

    cp = cps[0]

    assert cp.name == "Chn_ZERO_LENGTH"
    assert cp.bed_slope == 0.1234500000
    assert len(cp.survey_points) == 8
    assert cp.survey_points[0].root == (0, 1444.5898)
    assert cp.survey_points[1].root == (16.0000, 1440.5898)
    assert len(cp.roughness_zones) == 3
    assert cp.roughness_zones[0].root == (0, 0.12345000)

    assert str(cp) == dedent(
        """
    :ChannelProfile Chn_ZERO_LENGTH
      :Bedslope             0.12345

      :SurveyPoints
        0.0 1444.5898
        16.0 1440.5898
        16.2469 1440.5898
        16.2778 1440.4663
        16.3395 1440.4663
        16.3703 1440.5898
        16.6172 1440.5898
        32.6173 1444.5898
      :EndSurveyPoints

      :RoughnessZones
        0.0 0.12345
        16.2469 0.12345
        16.37035 0.12345
      :EndRoughnessZones
    :EndChannelProfile
    """
    )
