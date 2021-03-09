import datetime as dt
import os
import tempfile
import zipfile
from dataclasses import replace

import numpy as np
import pytest
import xarray as xr

from ravenpy.models import (
    GR4JCN,
    GR4JCN_OST,
    HBVEC,
    HBVEC_OST,
    HMETS,
    HMETS_OST,
    MOHYSE,
    MOHYSE_OST,
    Raven,
    Routing,
)
from ravenpy.models.commands import (
    GriddedForcingCommand,
    HRUStateVariableTableCommand,
    LandUseClassesCommand,
    ObservationDataCommand,
    SBGroupPropertyMultiplierCommand,
    SoilClassesCommand,
    SoilProfilesCommand,
    VegetationClassesCommand,
)
from ravenpy.utilities.testdata import get_local_testdata

from .common import _convert_2d

TS = get_local_testdata(
    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
)


@pytest.fixture
def input2d(tmpdir):
    """Convert 1D input to 2D output by copying all the time series along a new region dimension."""
    ds = _convert_2d(TS)
    fn_out = os.path.join(tmpdir, "input2d.nc")
    ds.to_netcdf(fn_out)
    return fn_out


def test_race():
    model1 = GR4JCN()
    model1.rvi.suppress_output = True
    model2 = GR4JCN()
    ost = GR4JCN_OST()

    assert model1.rvi.suppress_output.startswith(":SuppressOutput")
    assert model2.rvi.suppress_output == ""
    assert ost.rvi.suppress_output.startswith(":SuppressOutput")


class TestGR4JCN:
    def test_simple(self):
        model = GR4JCN(tempfile.mkdtemp())

        model.rvi.start_date = dt.datetime(2000, 1, 1)
        model.rvi.end_date = dt.datetime(2002, 1, 1)
        model.rvi.run_name = "test"

        model.rvh.name = "Salmon"
        model.rvh.area = "4250.6"
        model.rvh.elevation = "843.0"
        model.rvh.latitude = 54.4848
        model.rvh.longitude = -123.3659

        model.rvt.pr.deaccumulate = False

        model.rvp.params = model.params(0.529, -3.396, 407.29, 1.072, 16.9, 0.947)
        assert model.rvi.suppress_output == ""
        model(TS)

        d = model.diagnostics
        # yields NSE=0.???? for full period 1954-2010
        assert model.rvi.calendar == "GREGORIAN"
        # Check parser
        assert 1 in model.solution["HRUStateVariableTable"]["data"]

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -0.117301, 2)

        hds = model.q_sim
        assert hds.attrs["long_name"] == "Simulated outflows"

        # Check attributes
        assert model.hydrograph.attrs["model_id"] == "gr4jcn"

    def test_tags(self):
        model = GR4JCN(tempfile.mkdtemp())

        tags = model.tags
        assert "run_name" in tags

    def test_rvobjs(self):
        model = GR4JCN(tempfile.mkdtemp())
        a = model.rvobjs
        assert a

    def test_assign(self):
        model = GR4JCN()
        model.assign("run_name", "test")
        assert model.rvi.run_name == "test"

        model.assign("params", np.array([0.529, -3.396, 407.29, 1.072, 16.9, 0.947]))
        assert model.rvp.params.GR4J_X1 == 0.529

        model.assign("params", [0.529, -3.396, 407.29, 1.072, 16.9, 0.947])
        assert model.rvp.params.GR4J_X1 == 0.529

        model.assign("params", (0.529, -3.396, 407.29, 1.072, 16.9, 0.947))
        assert model.rvp.params.GR4J_X1 == 0.529

    def test_run(self):
        model = GR4JCN()
        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
            suppress_output=False,
        )
        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -0.117301, 2)

    # @pytest.mark.skip
    def test_overwrite(self):
        model = GR4JCN()
        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
        )
        assert model.rvi.suppress_output == ""

        qsim1 = model.q_sim.copy(deep=True)
        m1 = qsim1.mean()

        # This is only needed temporarily while we fix this: https://github.com/CSHS-CWRA/RavenPy/issues/4
        # Please remove when fixed!
        model.hydrograph.close()  # Needed with xarray 0.16.1

        model(TS, params=(0.5289, -3.397, 407.3, 1.071, 16.89, 0.948), overwrite=True)

        qsim2 = model.q_sim.copy(deep=True)
        m2 = qsim2.mean()

        # This is only needed temporarily while we fix this: https://github.com/CSHS-CWRA/RavenPy/issues/4
        # Please remove when fixed!
        model.hydrograph.close()  # Needed with xarray 0.16.1

        assert m1 != m2

        np.testing.assert_almost_equal(m1, m2, 1)

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -0.117315, 2)

        # Set initial conditions explicitly
        model(
            TS,
            end_date=dt.datetime(2001, 2, 1),
            hru_state=HRUStateVariableTableCommand.Record(soil0=0),
            overwrite=True,
        )
        assert model.q_sim.isel(time=1).values[0] < qsim2.isel(time=1).values[0]

    def test_resume(self):
        model_ab = GR4JCN()
        kwargs = dict(
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
        )
        # Reference run
        model_ab(
            TS,
            run_name="run_ab",
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2001, 1, 1),
            **kwargs,
        )

        model_a = GR4JCN()
        model_a(
            TS,
            run_name="run_a",
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2000, 7, 1),
            **kwargs,
        )

        # Path to solution file from run A
        rvc = model_a.outputs["solution"]

        # Resume with final state from live model
        model_a.resume()

        model_a(
            TS,
            run_name="run_2",
            start_date=dt.datetime(2000, 7, 1),
            end_date=dt.datetime(2001, 1, 1),
            **kwargs,
        )

        for key in ["Soil Water[0]", "Soil Water[1]"]:
            np.testing.assert_array_almost_equal(
                model_a.storage[key] - model_ab.storage[key], 0, 5
            )

        # Resume with final state from saved solution file
        model_b = GR4JCN()
        model_b.resume(
            rvc
        )  # <--------- And this is how you feed it to a brand new model.
        model_b(
            TS,
            run_name="run_2",
            start_date=dt.datetime(2000, 7, 1),
            end_date=dt.datetime(2001, 1, 1),
            **kwargs,
        )

        for key in ["Soil Water[0]", "Soil Water[1]"]:
            np.testing.assert_array_almost_equal(
                model_b.storage[key] - model_ab.storage[key], 0, 5
            )

        # model.solution loads the solution in a dictionary. I expected the variables to be identical,
        # but some atmosphere related attributes are way off. Is it possible that `ATMOSPHERE` and `ATMOS_PRECIP` are
        # cumulative sums of precipitation over the run ?
        # assert model_b.solution == model_ab.solution # This does not work. Atmosphere attributes are off.

    def test_resume_earlier(self):
        """Check that we can resume a run with the start date set at another date than the time stamp in the
        solution."""
        kwargs = dict(
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
        )
        # Reference run
        model = GR4JCN()
        model(
            TS,
            run_name="run_a",
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2000, 2, 1),
            **kwargs,
        )

        s_a = model.storage["Soil Water[0]"].isel(time=-1)

        # Path to solution file from run A
        rvc = model.outputs["solution"]

        # Resume with final state from live model
        # We have two options to do this:
        # 1. Replace model template by solution file as is: model.resume()
        # 2. Replace variable in RVC class by parsed values: model.rvc.parse(rvc.read_text())
        # I think in many cases option 2 will prove simpler.

        model.rvc.parse(rvc.read_text())

        model(
            TS,
            run_name="run_b",
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2000, 2, 1),
            **kwargs,
        )

        s_b = model.storage["Soil Water[0]"].isel(time=-1)
        assert s_a != s_b

    def test_update_soil_water(self):
        kwargs = dict(
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
        )
        # Reference run
        model = GR4JCN()
        model(
            TS,
            run_name="run_a",
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2000, 2, 1),
            **kwargs,
        )

        s_0 = float(model.storage["Soil Water[0]"].isel(time=-1).values)
        s_1 = float(model.storage["Soil Water[1]"].isel(time=-1).values)

        hru_state = replace(model.rvc.hru_state, soil0=s_0, soil1=s_1)

        model(
            TS,
            run_name="run_b",
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2000, 2, 1),
            hru_state=hru_state,
            **kwargs,
        )

        assert s_0 != model.storage["Soil Water[0]"].isel(time=-1)
        assert s_1 != model.storage["Soil Water[1]"].isel(time=-1)

    def test_version(self):
        model = Raven()
        assert model.version == "3.0.1"

        model = GR4JCN()
        assert model.version == "3.0.1"

    def test_parallel_params(self):
        model = GR4JCN()
        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=[
                (0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
                (0.528, -3.4, 407.3, 1.07, 17, 0.95),
            ],
            suppress_output=False,
        )

        assert len(model.diagnostics) == 2
        assert model.hydrograph.dims["params"] == 2
        z = zipfile.ZipFile(model.outputs["rv_config"])
        assert len(z.filelist) == 10

    def test_parallel_basins(self, input2d):
        ts = input2d
        model = GR4JCN()
        model(
            ts,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=[0.529, -3.396, 407.29, 1.072, 16.9, 0.947],
            nc_index=[0, 0],
            name=["basin1", "basin2"],
            suppress_output=False,
        )

        assert len(model.diagnostics) == 2
        assert len(model.hydrograph.nbasins) == 2
        np.testing.assert_array_equal(
            model.hydrograph.basin_name[:], ["basin1", "basin2"]
        )
        z = zipfile.ZipFile(model.outputs["rv_config"])
        assert len(z.filelist) == 10


class TestGR4JCN_OST:
    def test_simple(self):
        model = GR4JCN_OST()
        params = (0.529, -3.396, 407.29, 1.072, 16.9, 0.053)
        low = (0.01, -15.0, 10.0, 0.0, 1.0, 0.0)
        high = (2.5, 10.0, 700.0, 7.0, 30.0, 1.0)

        model(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            lowerBounds=low,
            upperBounds=high,
            algorithm="DDS",
            random_seed=0,
            max_iterations=10,
        )

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], 0.50717, 4)

        # Random number seed: 123
        # Budget:             10
        # Algorithm:          DDS
        # :StartDate          1954-01-01 00:00:00
        # :Duration           208
        opt_para = model.calibrated_params
        opt_func = model.obj_func

        np.testing.assert_almost_equal(
            opt_para,
            [2.424726, 3.758972, 204.3856, 5.866946, 16.60408, 0.3728098],
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            -0.50717,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        # # Random number seed: 123
        # # Budget:             50
        # # Algorithm:          DDS
        # # :StartDate          1954-01-01 00:00:00
        # # :Duration           20819
        # np.testing.assert_almost_equal( opt_para, [0.3243268,3.034247,407.2890,2.722774,12.18124,0.9468769], 4,
        #                                 err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal( opt_func, -0.5779910, 4,
        #                                 err_msg='calibrated NSE is not matching expected value')
        gr4j = GR4JCN()
        gr4j(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=model.calibrated_params,
        )
        np.testing.assert_almost_equal(
            gr4j.diagnostics["DIAG_NASH_SUTCLIFFE"], d["DIAG_NASH_SUTCLIFFE"]
        )


class TestHMETS:
    def test_simple(self):
        model = HMETS()
        params = (
            9.5019,
            0.2774,
            6.3942,
            0.6884,
            1.2875,
            5.4134,
            2.3641,
            0.0973,
            0.0464,
            0.1998,
            0.0222,
            -1.0919,
            2.6851,
            0.3740,
            1.0000,
            0.4739,
            0.0114,
            0.0243,
            0.0069,
            310.7211,
            916.1947,
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            suppress_output=True,
        )

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -3.0132, 4)


class TestHMETS_OST:
    def test_simple(self):
        model = HMETS_OST()
        params = (
            9.5019,
            0.2774,
            6.3942,
            0.6884,
            1.2875,
            5.4134,
            2.3641,
            0.0973,
            0.0464,
            0.1998,
            0.0222,
            -1.0919,
            2.6851,
            0.3740,
            1.0000,
            0.4739,
            0.0114,
            0.0243,
            0.0069,
            310.7211,
            916.1947,
        )
        low = (
            0.3,
            0.01,
            0.5,
            0.15,
            0.0,
            0.0,
            -2.0,
            0.01,
            0.0,
            0.01,
            0.005,
            -5.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.00001,
            0.0,
            0.00001,
            0.0,
            0.0,
        )
        high = (
            20.0,
            5.0,
            13.0,
            1.5,
            20.0,
            20.0,
            3.0,
            0.2,
            0.1,
            0.3,
            0.1,
            2.0,
            5.0,
            1.0,
            3.0,
            1.0,
            0.02,
            0.1,
            0.01,
            0.5,
            2.0,
        )

        model(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            lowerBounds=low,
            upperBounds=high,
            algorithm="DDS",
            random_seed=0,
            max_iterations=10,
        )

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -1.43474, 4)

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        # # Random number seed: 123
        # # Budget:             50
        # # Algorithm:          DDS
        # # :StartDate          1954-01-01 00:00:00
        # # :Duration           20819
        # np.testing.assert_almost_equal( opt_para, [0.3243268,3.034247,407.2890,2.722774,12.18124,0.9468769], 4,
        #                                 err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal( opt_func, -0.5779910, 4,
        #                                 err_msg='calibrated NSE is not matching expected value')
        #
        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #

        expected_value = [
            1.777842e01,
            3.317211e00,
            5.727342e00,
            1.419491e00,
            1.382141e01,
            1.637954e01,
            7.166296e-01,
            1.389346e-01,
            2.620464e-02,
            2.245525e-01,
            2.839426e-02,
            -2.003810e00,
            9.479623e-01,
            4.803857e-01,
            2.524914e00,
            4.117232e-01,
            1.950058e-02,
            4.494123e-02,
            1.405815e-03,
            2.815803e-02,
            1.007823e00,
        ]
        np.testing.assert_almost_equal(
            opt_para,
            expected_value,
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            1.43474,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-hmets
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [5.008045E+00, 7.960246E-02, 4.332698E+00, 4.978125E-01,
        #                                           1.997029E+00, 6.269773E-01, 1.516961E+00, 8.180383E-02,
        #                                           6.730663E-02, 2.137822E-02, 2.097163E-02, 1.773348E+00,
        #                                           3.036039E-01, 1.928524E-02, 1.758471E+00, 8.942299E-01,
        #                                           8.741980E-03, 5.036474E-02, 9.465804E-03, 1.851839E-01,
        #                                           1.653934E-01, 2.624006E+00, 8.868485E-02, 9.259195E+01,
        #                                           8.269670E+01], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -6.350490E-01, 4,
        #                                err_msg='calibrated NSE is not matching expected value')
        hmets = HMETS()
        hmets(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=model.calibrated_params,
        )

        np.testing.assert_almost_equal(
            hmets.diagnostics["DIAG_NASH_SUTCLIFFE"], d["DIAG_NASH_SUTCLIFFE"], 4
        )


class TestMOHYSE:
    def test_simple(self):
        model = MOHYSE()
        params = (
            1.0,
            0.0468,
            4.2952,
            2.658,
            0.4038,
            0.0621,
            0.0273,
            0.0453,
            0.9039,
            5.6167,
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            suppress_output=True,
        )

        d = model.diagnostics
        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], 0.194612, 4)


class TestMOHYSE_OST:
    def test_simple(self):
        model = MOHYSE_OST()
        params = (
            1.0,
            0.0468,
            4.2952,
            2.658,
            0.4038,
            0.0621,
            0.0273,
            0.0453,
            0.9039,
            5.6167,
        )

        low_p = (0.01, 0.01, 0.01, -5.00, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
        high_p = (20.0, 1.0, 20.0, 5.0, 0.5, 1.0, 1.0, 1.0, 15.0, 15.0)

        model(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            lowerBounds=low_p,
            upperBounds=high_p,
            algorithm="DDS",
            random_seed=0,
            max_iterations=10,
        )

        d = model.diagnostics
        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], 0.3826810, 4)

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        # # Random number seed: 123
        # # Budget:             50
        # # Algorithm:          DDS
        # # :StartDate          1954-01-01 00:00:00
        # # :Duration           20819
        # np.testing.assert_almost_equal( opt_para, [0.3243268,3.034247,407.2890,2.722774,12.18124,0.9468769], 4,
        #                                 err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal( opt_func, -0.5779910, 4,
        #                                 err_msg='calibrated NSE is not matching expected value')
        #
        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #
        np.testing.assert_almost_equal(
            opt_para,
            [
                7.721801e00,
                8.551484e-01,
                1.774571e01,
                1.627677e00,
                7.702450e-02,
                9.409600e-01,
                6.941596e-01,
                8.207870e-01,
                8.154455e00,
                1.018226e01,
            ],
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            -0.3826810,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-mohyse
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [1.517286E+01, 7.112556E-01, 1.981243E+01, -4.193046E+00,
        #                                           1.791486E-01, 9.774897E-01, 5.353541E-01, 6.686806E-01,
        #                                           1.040908E+01, 1.132304E+01, 8.831552E-02], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -0.3857010, 4,
        #                                err_msg='calibrated NSE is not matching expected value')
        mohyse = MOHYSE()
        mohyse(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=model.calibrated_params,
        )

        np.testing.assert_almost_equal(
            mohyse.diagnostics["DIAG_NASH_SUTCLIFFE"], d["DIAG_NASH_SUTCLIFFE"], 4
        )


class TestHBVEC:
    def test_simple(self):
        model = HBVEC()
        params = (
            0.05984519,
            4.072232,
            2.001574,
            0.03473693,
            0.09985144,
            0.506052,
            3.438486,
            38.32455,
            0.4606565,
            0.06303738,
            2.277781,
            4.873686,
            0.5718813,
            0.04505643,
            0.877607,
            18.94145,
            2.036937,
            0.4452843,
            0.6771759,
            1.141608,
            1.024278,
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            suppress_output=True,
        )

        d = model.diagnostics
        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], 0.0186633, 4)

    def test_evap(self):
        model = HBVEC()
        params = (
            0.05984519,
            4.072232,
            2.001574,
            0.03473693,
            0.09985144,
            0.506052,
            3.438486,
            38.32455,
            0.4606565,
            0.06303738,
            2.277781,
            4.873686,
            0.5718813,
            0.04505643,
            0.877607,
            18.94145,
            2.036937,
            0.4452843,
            0.6771759,
            1.141608,
            1.024278,
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            suppress_output=True,
            evaporation="PET_OUDIN",
            ow_evaporation="PET_OUDIN",
        )


class TestHBVEC_OST:
    def test_simple(self):
        model = HBVEC_OST()
        params = (
            0.05984519,
            4.072232,
            2.001574,
            0.03473693,
            0.09985144,
            0.506052,
            3.438486,
            38.32455,
            0.4606565,
            0.06303738,
            2.277781,
            4.873686,
            0.5718813,
            0.04505643,
            0.877607,
            18.94145,
            2.036937,
            0.4452843,
            0.6771759,
            1.141608,
            1.024278,
        )

        low = (
            -3.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.3,
            0.0,
            0.0,
            0.01,
            0.05,
            0.01,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.01,
            0.0,
            0.05,
            0.8,
            0.8,
        )
        high = (
            3.0,
            8.0,
            8.0,
            0.1,
            1.0,
            1.0,
            7.0,
            100.0,
            1.0,
            0.1,
            6.0,
            5.0,
            5.0,
            0.2,
            1.0,
            30.0,
            3.0,
            2.0,
            1.0,
            1.5,
            1.5,
        )

        model(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=params,
            lowerBounds=low,
            upperBounds=high,
            algorithm="DDS",
            random_seed=0,
            max_iterations=10,
        )

        d = model.diagnostics
        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -2.25991e-01, 4)

        opt_para = model.calibrated_params
        opt_func = model.obj_func

        # Random number seed: 123                         #
        # Budget:             10                          #      This is the setup used for testing:
        # Algorithm:          DDS                         #         shorter sim-period and lower budget
        # :StartDate          1954-01-01 00:00:00         #      First tested that example below matches
        # :Duration           208                         #
        np.testing.assert_almost_equal(
            opt_para,
            [
                -8.317931e-01,
                4.072232e00,
                2.001574e00,
                5.736299e-03,
                9.985144e-02,
                4.422529e-01,
                3.438486e00,
                8.055843e01,
                4.440133e-01,
                8.451082e-02,
                2.814201e00,
                7.327970e-01,
                1.119773e00,
                1.161223e-03,
                4.597179e-01,
                1.545857e01,
                1.223865e00,
                4.452843e-01,
                9.492006e-01,
                9.948123e-01,
                1.110682e00,
            ],
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            2.25991e-01,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        # # Random number seed: 123                       #
        # # Budget:             50                        #      This is the setup in the Wiki:
        # # Algorithm:          DDS                       #      https://github.com/Ouranosinc/raven/wiki/
        # # :StartDate          1954-01-01 00:00:00       #      Technical-Notes#example-setups-for-environment-
        # # :Duration           20819                     #
        # np.testing.assert_almost_equal(opt_para, [5.984519E-02, 4.072232E+00, 2.001574E+00, 3.473693E-02,
        #                                           9.985144E-02, 5.060520E-01, 2.944343E+00, 3.832455E+01,
        #                                           4.606565E-01, 6.303738E-02, 2.277781E+00, 4.873686E+00,
        #                                           5.718813E-01, 4.505643E-02, 8.776511E-01, 1.894145E+01,
        #                                           2.036937E+00, 4.452843E-01, 6.771759E-01, 1.206053E+00,
        #                                           1.024278E+00], 4,
        #                                err_msg='calibrated parameter set is not matching expected value')
        # np.testing.assert_almost_equal(opt_func, -6.034670E-01, 4,
        #                                err_msg='calibrated NSE is not matching expected value')
        hbvec = HBVEC()
        hbvec(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            area=4250.6,
            elevation=843.0,
            latitude=54.4848,
            longitude=-123.3659,
            params=model.calibrated_params,
        )

        np.testing.assert_almost_equal(
            hbvec.diagnostics["DIAG_NASH_SUTCLIFFE"], d["DIAG_NASH_SUTCLIFFE"], 4
        )


class TestRouting:
    importers = pytest.importorskip("ravenpy.models.importers")

    def test_lievre_tutorial(self):
        """
        This test reproduces the Lievre tutorial setup:

        http://raven.uwaterloo.ca/files/RavenTutorial6.zip

        """

        ###############
        # Input files #
        ###############

        routing_product_shp_path = get_local_testdata(
            "raven-routing-sample/finalcat_hru_info.zip"
        )

        vic_streaminputs_nc_path = get_local_testdata(
            "raven-routing-sample/VIC_streaminputs.nc"
        )
        vic_temperatures_nc_path = get_local_testdata(
            "raven-routing-sample/VIC_temperatures.nc"
        )

        observation_data_nc_path = get_local_testdata(
            "raven-routing-sample/WSC02LE024.nc"
        )

        #########
        # Model #
        #########

        model = Routing("/tmp/ravenpy_routing_emu_dev/")

        #######
        # RVI #
        #######

        streaminputs = xr.open_dataset(vic_streaminputs_nc_path)

        start = streaminputs.indexes["time"][0]
        end = streaminputs.indexes["time"][-4]  # to match the tutorial end date

        model.rvi.start_date = start.to_pydatetime()
        model.rvi.end_date = end.to_pydatetime()

        # Raven will use 24h even though the NC inputs are 6h
        model.rvi.time_step = "24:00:00"

        model.rvi.evaluation_metrics = "NASH_SUTCLIFFE PCT_BIAS KLING_GUPTA"

        #######
        # RVH #
        #######

        rvh_importer = self.importers.RoutingProductShapefileImporter(
            routing_product_shp_path, hru_aspect_convention="ArcGIS"
        )
        rvh_config = rvh_importer.extract()
        channel_profiles = rvh_config.pop("channel_profiles")

        gauge = [sb for sb in rvh_config["subbasins"] if sb.gauged]
        assert len(gauge) == 1
        gauge = gauge[0]

        model.rvh.update(rvh_config)

        model.rvh["sb_group_property_multipliers"] = (
            SBGroupPropertyMultiplierCommand("Land", "MANNINGS_N", 1.0),
            SBGroupPropertyMultiplierCommand("Lakes", "RESERVOIR_CREST_WIDTH", 1.0),
        )

        #######
        # RVP #
        #######

        # The labels used for the following commands ("Lake_Soil_Lake_HRU", "Soil_Land_HRU", etc)
        # must correspond to the values of certain fields of the Routing Product:
        # LAND_USE_C, VEG_C, SOIL_PROF

        model.rvp.avg_annual_runoff = 594
        model.rvp.soil_classes = [SoilClassesCommand.Record("AQUIFER")]
        model.rvp.soil_profiles = [
            SoilProfilesCommand.Record("Lake_Soil_Lake_HRU", 1, "AQUIFER", 5),
            SoilProfilesCommand.Record("Soil_Land_HRU", 1, "AQUIFER", 5),
        ]
        model.rvp.vegetation_classes = [
            VegetationClassesCommand.Record("Veg_Land_HRU", 25, 5.0, 5.0),
            VegetationClassesCommand.Record("Veg_Lake_HRU", 0, 0, 0),
        ]
        model.rvp.land_use_classes = [
            LandUseClassesCommand.Record("Landuse_Land_HRU", 0, 1),
            LandUseClassesCommand.Record("Landuse_Lake_HRU", 0, 0),
        ]
        model.rvp.channel_profiles = channel_profiles

        #######
        # RVT #
        #######

        streaminputs_importer = self.importers.RoutingProductGridWeightImporter(
            vic_streaminputs_nc_path, routing_product_shp_path
        )

        temperatures_importer = self.importers.RoutingProductGridWeightImporter(
            vic_temperatures_nc_path, routing_product_shp_path
        )

        streaminputs_gf = GriddedForcingCommand(
            name="StreamInputs",
            forcing_type="PRECIP",
            file_name_nc=vic_streaminputs_nc_path.name,
            var_name_nc="Streaminputs",
            dim_names_nc=("lon_dim", "lat_dim", "time"),
            linear_transform=(4.0, 0),
            grid_weights=streaminputs_importer.extract(),
        )

        temperatures_gf = GriddedForcingCommand(
            name="AverageTemp",
            forcing_type="TEMP_AVE",
            file_name_nc=vic_temperatures_nc_path.name,
            var_name_nc="Avg_temp",
            dim_names_nc=("lon_dim", "lat_dim", "time"),
            grid_weights=temperatures_importer.extract(),
        )

        model.rvt.gridded_forcings = [streaminputs_gf, temperatures_gf]

        model.rvt.observation_data = ObservationDataCommand(
            data_type="HYDROGRAPH",
            subbasin_id=gauge.subbasin_id,
            units="m3/s",
            file_name_nc=observation_data_nc_path.name,
            var_name_nc="Q",
            dim_names_nc=("nstations", "time"),
            station_idx=1,
        )

        #############
        # Run model #
        #############

        model(
            [
                vic_streaminputs_nc_path,
                vic_temperatures_nc_path,
                observation_data_nc_path,
            ],
            overwrite=True,
        )

        ##########
        # Verify #
        ##########

        assert len(model.hydrograph.time) == (end - start).days + 1

        assert model.hydrograph.basin_name.item() == gauge.name

        csv_lines = model.outputs["diagnostics"].read_text().split("\n")
        assert csv_lines[1].split(",")[:-1] == [
            "HYDROGRAPH_ALL",
            observation_data_nc_path.name,
            "0.311049",
            "-11.9514",
            "0.471256",
        ]

        for d, q_sim in [
            (0, 85.97869699520872),
            (1000, 78.9516829670743),
            (2000, 70.29412760595676),
            (3000, 44.711237489482755),
            (4000, 129.98874279175033),
        ]:
            assert model.hydrograph.q_sim[d].item() == pytest.approx(q_sim)
