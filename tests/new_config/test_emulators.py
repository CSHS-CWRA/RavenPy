import datetime as dt

import numpy as np
import pytest
import xarray as xr

from ravenpy import Emulator
from ravenpy.new_config import commands as rc
from ravenpy.new_config.emulators import GR4JCN
from ravenpy.new_config.rvs import Config
from ravenpy.new_config.utils import get_average_annual_runoff

# Expected NSE for emulator configuration from the `config_rv` test fixture.
NSE = {
    "GR4JCN": -0.117301,
    "HMETS": -3.0132,
    "Mohyse": 0.194612,
    "HBVEC": 0.0186633,
    "CanadianShield": 0.39602,
    "HYPR": 0.685188,
    "SACSMA": -0.0382907,
    "Blended": -0.913785,
}


def test_rv(config_rv):
    """Test the model configuration can be written to disk."""
    name, path = config_rv
    assert (path / "test.rvi").exists()


def test_run(numeric_config, tmp_path):
    """Test that the emulator actually runs and returns the expected NSE."""
    name, conf = numeric_config
    assert conf.__config__.allow_mutation

    e = Emulator(config=conf, workdir=tmp_path)
    e.write_rv()
    out = e.run()

    d = out.diagnostics
    np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], NSE[name], 4)

    if name == "CanadianShield":
        pytest.skip("Missing solution due to SuppressOutput.")

    # Start new simulation with final state from initial run.
    new = e.resume()
    new.start_date = conf.end_date
    new.end_date = dt.datetime(2002, 1, 7)
    new.run_name = None

    e2 = Emulator(config=new, workdir=tmp_path / "resumed")
    out = e2.run()
    assert e2.modelname == "raven"
    assert isinstance(out.hydrograph, xr.Dataset)
    assert isinstance(out.storage, xr.Dataset)


@pytest.mark.skip("Need to find a clean way to freeze emulator config instance.")
def test_emulator_config_is_read_only(dummy_config, tmp_path):
    cls, P = dummy_config

    e = Emulator(config=cls(), workdir=tmp_path)

    # The emulator configuration should be read-only.
    with pytest.raises(TypeError):
        e.config.run_name = "Renamed"


def test_no_run_name(dummy_config, tmp_path):
    cls, P = dummy_config
    conf = cls(params=[0.5])

    paths = conf.write_rv(tmp_path, modelname="ham")
    assert "ham.rvi" in str(paths["rvi"])

    paths = conf.write_rv(tmp_path, overwrite=True)
    assert "raven.rvi" in str(paths["rvi"])


@pytest.mark.slow
@pytest.mark.online
def test_run_with_dap_link(minimal_emulator, tmp_path):
    """Test Raven with DAP link instead of local netCDF file."""
    # Link to THREDDS Data Server netCDF testdata
    TDS = "https://pavics.ouranos.ca/twitcher/ows/proxy/thredds/dodsC/birdhouse/testdata/raven"
    fn = f"{TDS}/raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"

    alt_names = {
        "RAINFALL": "rain",
        "TEMP_MIN": "tmin",
        "TEMP_MAX": "tmax",
        "PET": "pet",
        "HYDROGRAPH": "qobs",
        "SNOWFALL": "snow",
    }

    conf = minimal_emulator
    conf.gauge = [rc.Gauge.from_nc(fn, alt_names=alt_names)]

    out = Emulator(conf, workdir=tmp_path).run()


def test_routing(get_local_testdata):
    """We need at least 2 subbasins to activate routing."""
    from ravenpy.new_config.emulators.gr4jcn import P

    # Salmon catchment is now split into land- and lake-part.
    # The areas do not sum up to overall area of 4250.6 [km2].
    # This is the reason the "test_routing" will give different
    # results compared to "test_run". The "salmon_land_hru"
    # however is kept at the overall area of 4250.6 [km2] such
    # that other tests still obtain same results as before.

    salmon_land_hru_1 = dict(
        area=4250.6,
        elevation=843.0,
        latitude=54.4848,
        longitude=-123.3659,
        hru_type="land",
    )
    salmon_lake_hru_1 = dict(
        area=100.0, elevation=839.0, latitude=54.0, longitude=-123.4, hru_type="lake"
    )
    salmon_land_hru_2 = dict(
        area=2000.0,
        elevation=835.0,
        latitude=54.123,
        longitude=-123.4234,
        hru_type="land",
    )

    salmon_river = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    ts_2d = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_2d.nc"
    )

    #########
    # R V P #
    #########

    parameters = [0.529, -3.396, 407.29, 1.072, 16.9, 0.947]

    #########
    # R V H #
    #########

    # Here we assume that we have two subbasins. The first one (subbasin_id=10)
    # has a lake (hru_id=2; area-100km2) and the rest is covered by land (hru_id=1;
    # area=4250.6km2). The second subbasin (subbasin_id=20) does not contain a
    # lake and is hence only land (hru_id=3; area=2000km2).
    #
    # Later the routing product will tell us which basin flows into which. Here
    # we assume that the first subbasin (subbasin_id=10) drains into the second
    # (subbasin_id=20). At the outlet of this second one we have an observation
    # station (see :ObservationData in RVT). We will compare these observations
    # with the simulated streamflow. That is the reason why "gauged=True" for
    # the second basin.

    # HRU IDs are 1 to 3
    # Sub-basin IDs are 10 and 20 (not 1 and 2), to help disambiguate
    hru1 = dict(hru_id=1, subbasin_id=10, **salmon_land_hru_1)
    hru2 = dict(hru_id=2, subbasin_id=10, **salmon_lake_hru_1)
    hru3 = dict(hru_id=3, subbasin_id=20, **salmon_land_hru_2)
    hrus = [hru1, hru2, hru3]

    hs = [
        rc.HRUState(hru_id=1, data={"SOIL[0]": P.GR4J_X1 * 1000 / 2, "SOIL[1]": 15}),
        rc.HRUState(hru_id=2, data={"SOIL[0]": 0, "SOIL[1]": 0}),
        rc.HRUState(hru_id=3, data={"SOIL[0]": P.GR4J_X1 * 1000 / 2, "SOIL[1]": 15}),
    ]

    # gauged = False:
    # Usually this output would only be written for user's convenience.
    # There is usually no observation of streamflow available within
    # catchments; only at the outlet. That's most commonly the reason
    # why a catchment is defined as it is defined.
    sub1 = dict(
        name="upstream",
        subbasin_id=10,
        downstream_id=20,
        profile="chn_10",
        gauged=False,
    )
    # gauged = True:
    # Since this is the outlet, this would usually be what we calibrate
    # against (i.e. we try to match this to Qobs).
    sub2 = dict(
        name="downstream",
        subbasin_id=20,
        downstream_id=-1,
        profile="chn_20",
        gauged=True,
    )

    # Set basin state variables to 0
    bs = [rc.BasinIndex(sb_id=10), rc.BasinIndex(sb_id=20)]

    SB_groups = [
        rc.SubBasinGroup(name="Land", sb_ids=[10]),
        rc.SubBasinGroup(name="Lakes", sb_ids=[20]),
    ]
    SBP = [
        rc.SBGroupPropertyMultiplier(
            group_name="Land", parameter_name="MANNINGS_N", mult=1.0
        ),
        rc.SBGroupPropertyMultiplier(
            group_name="Lakes", parameter_name="RESERVOIR_CREST_WIDTH", mult=1.0
        ),
    ]

    total_area_in_m2 = sum(hru["area"] for hru in hrus) * 1000 * 1000
    avg_annual_runoff = get_average_annual_runoff(ts_2d, total_area_in_m2)
    np.testing.assert_almost_equal(avg_annual_runoff, 139.5407534171111)

    r"""
    # These channel profiles describe the geometry of the actual river crossection.
    # The eight points (x) to describe the following geometry are given in each
    # profile:
    #
    # ----x                                     x---
    #      \           FLOODPLAIN             /
    #       x----x                     x----x
    #             \                  /
    #               \   RIVERBED   /
    #                 x----------x
    #
    """
    channel_profiles = [
        rc.ChannelProfile(
            name="chn_10",
            bed_slope=7.62066e-05,
            survey_points=[
                (0, 463.647),
                (16.0, 459.647),
                (90.9828, 459.647),
                (92.9828, 458.647),
                (126.4742, 458.647),
                (128.4742, 459.647),
                (203.457, 459.647),
                (219.457, 463.647),
            ],
            roughness_zones=[
                (0, 0.0909167),
                (90.9828, 0.035),
                (128.4742, 0.0909167),
            ],
        ),
        rc.ChannelProfile(
            name="chn_20",
            bed_slope=9.95895e-05,
            survey_points=[
                (0, 450.657),
                (16.0, 446.657),
                (85.0166, 446.657),
                (87.0166, 445.657),
                (117.5249, 445.657),
                (119.5249, 446.657),
                (188.54149999999998, 446.657),
                (204.54149999999998, 450.657),
            ],
            roughness_zones=[
                (0, 0.0915769),
                (85.0166, 0.035),
                (119.5249, 0.0915769),
            ],
        ),
    ]
    #########
    # R V T #
    #########

    alt_names = {
        "RAINFALL": "rain",
        "TEMP_MIN": "tmin",
        "TEMP_MAX": "tmax",
        "SNOWFALL": "snow",
        "PET": "pet",
    }

    data_type = ["RAINFALL", "TEMP_MIN", "TEMP_MAX", "SNOWFALL", "PET"]
    gw = rc.GridWeights(
        number_hrus=3,
        number_grid_cells=1,
        # Here we have a special case: station is 0 for every row because the example NC
        # has only one region/station (which is column 0)
        data=((1, 0, 1.0), (2, 0, 1.0), (3, 0, 1.0)),
    )

    sf = [
        rc.StationForcing.from_nc(ts_2d, typ, alt_names=alt_names[typ], GridWeights=gw)
        for typ in data_type
    ]

    obs = [rc.ObservationData.from_nc(salmon_river, alt_names="qobs", uid=20)]

    ###################
    # Configure model #
    ###################

    model = GR4JCN(
        StartDate=dt.datetime(2000, 1, 1),
        EndDate=dt.datetime(2002, 1, 1),
        RunName="test_gr4jcn_routing",
        Routing="ROUTE_DIFFUSIVE_WAVE",
        params=parameters,
        HRUs=hrus,
        HRUStateVariableTable=hs,
        SubBasins=[sub1, sub2],
        SubBasinGroup=SB_groups,
        SBGroupPropertyMultiplier=SBP,
        BasinStateVariables=bs,
        ChannelProfile=channel_profiles,
        StationForcing=sf,
        ObservationData=obs,
        GlobalParameter={"AVG_ANNUAL_RUNOFF": avg_annual_runoff},
        WriteForcingFunctions=True,
        UniformInitialConditions=None,
        EvaluationMetrics=("NASH_SUTCLIFFE",),
    )

    #############
    # Run model #
    #############

    out = Emulator(model).run()

    ###########
    # Verify  #
    ###########

    hds = out.hydrograph.q_sim

    assert len(hds.nbasins) == 1  # number of "gauged" basins is 1

    # We only have one SB with gauged=True, so the output has a single column.
    # The number of time steps simulated between (2000, 1, 1) and
    # (2002, 1, 1) is 732.
    assert hds.shape == (732, 1)

    # Check simulated streamflow at first three timesteps and three simulated
    # timesteps in the middle of the simulation period.
    dates = (
        "2000-01-01",
        "2000-01-02",
        "2000-01-03",
        "2001-01-30",
        "2001-01-31",
        "2001-02-01",
    )

    target_q_sim = [0.0, 0.304073, 0.980807, 17.54049, 17.409493, 17.437954]

    for t in range(6):
        np.testing.assert_almost_equal(
            hds.sel(nbasins=0, time=dates[t]), target_q_sim[t], 4
        )

    # For lumped GR4J model we have 1 subbasin and 1 HRU as well as no routing, no
    # channel profiles, and the area of the entire basin is 4250.6 [km2]. Comparison
    # of simulated and observed streamflow at outlet yielded:
    # np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -0.116971, 4)
    #
    # This is now a different value due to:
    # - basin we have here is larger (4250.6 [km2] + 100 [km2] + 2000.0 [km2])
    # - we do routing: so water from subbasin 1 needs some time to arrive at the
    #   outlet of subbasin 2
    d = out.diagnostics
    np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -0.0141168, 4)

    assert len(list(out.path.glob("*ForcingFunctions.nc"))) == 1
