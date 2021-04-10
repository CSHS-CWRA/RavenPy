import datetime as dt
import tempfile

import numpy as np

from ravenpy.config.commands import LU
from ravenpy.models import BLENDED, BLENDED_OST
from ravenpy.utilities.testdata import get_local_testdata

from .common import _convert_2d

TS = get_local_testdata(
    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
)

hru = BLENDED.ForestHRU(
    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659, slope=0.01234
)

lu = LU("FOREST", impermeable_frac=0.0, forest_coverage=0.02345)


class TestBLENDED:
    def test_simple(self):
        model = BLENDED()
        params = (
            2.930702e-02,  # par_x01
            2.211166e00,  # par_x02
            2.166229e00,  # par_x03
            0.0002254976,  # feed 10**par_x04; not par_x4=-3.646858076
            2.173976e01,  # par_x05
            1.565091e00,  # par_x06
            6.211146e00,  # par_x07
            9.313578e-01,  # par_x08
            3.486263e-02,  # par_x09
            0.251835,  # feed par_x09+x10; not par_x10=0.21697237
            0.0002279250,  # feed 10**par_x11; not par_x11=-3.642208036
            1.214339e00,  # par_x12
            4.736668e-02,  # par_x13
            0.2070342,  # feed par_x13+x14; not par_x14=0.15966752
            7.806324e-02,  # par_x15
            -1.336429e00,  # par_x16
            2.189741e-01,  # par_x17
            3.845617e00,  # par_x18
            2.950022e-01,  # par_x19
            4.827523e-01,  # par_x20
            4.099820e00,  # par_x21
            1.283144e01,  # par_x22
            5.937894e-01,  # par_x23
            1.651588e00,  # par_x24
            1.705806,  # feed par_x24+x25; not par_x25=0.054218
            3.719308e-01,  # par_x26
            7.121015e-02,  # par_x27
            1.906440e-02,  # par_x28
            4.080660e-01,  # par_x29
            9.415693e-01,  # par_x30
            -1.856108e00,  # par_x31
            2.356995e00,  # par_x32
            1.0e00,  # par_x33
            1.0e00,  # par_x34
            7.510967e-03,  # par_x35
            5.321608e-01,  # par_r01
            2.891977e-02,  # par_r02
            9.605330e-01,  # par_r03
            6.128669e-01,  # par_r04
            9.558293e-01,  # par_r05
            1.008196e-01,  # par_r06
            9.275730e-02,  # par_r07
            7.469583e-01,  # par_r08
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            hrus=(hru,),
            land_use_classes=(lu,),
            params=params,
            suppress_output=True,
        )

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -0.913785, 4)


class TestBLENDED_OST:
    def test_simple(self):
        model = BLENDED_OST()
        params = (
            2.930702e-02,  # par_x01
            2.211166e00,  # par_x02
            2.166229e00,  # par_x03
            0.0002254976,  # feed 10**par_x04; not par_x4=-3.646858076
            2.173976e01,  # par_x05
            1.565091e00,  # par_x06
            6.211146e00,  # par_x07
            9.313578e-01,  # par_x08
            3.486263e-02,  # par_x09
            0.251835,  # feed par_x09+x10; not par_x10=0.21697237
            0.0002279250,  # feed 10**par_x11; not par_x11=-3.642208036
            1.214339e00,  # par_x12
            4.736668e-02,  # par_x13
            0.2070342,  # feed par_x13+x14; not par_x14=0.15966752
            7.806324e-02,  # par_x15
            -1.336429e00,  # par_x16
            2.189741e-01,  # par_x17
            3.845617e00,  # par_x18
            2.950022e-01,  # par_x19
            4.827523e-01,  # par_x20
            4.099820e00,  # par_x21
            1.283144e01,  # par_x22
            5.937894e-01,  # par_x23
            1.651588e00,  # par_x24
            1.705806,  # feed par_x24+x25; not par_x25=0.054218
            3.719308e-01,  # par_x26
            7.121015e-02,  # par_x27
            1.906440e-02,  # par_x28
            4.080660e-01,  # par_x29
            9.415693e-01,  # par_x30
            -1.856108e00,  # par_x31
            2.356995e00,  # par_x32
            1.110496e00,  # par_x33
            1.042556e00,  # par_x34
            7.510967e-03,  # par_x35
            5.321608e-01,  # par_r01
            2.891977e-02,  # par_r02
            9.605330e-01,  # par_r03
            6.128669e-01,  # par_r04
            9.558293e-01,  # par_r05
            1.008196e-01,  # par_r06
            9.275730e-02,  # par_r07
            7.469583e-01,  # par_r08
        )
        low = (
            0.0,  # par_x01
            0.1,  # par_x02
            0.5,  # par_x03
            -5.0,  # 10**par_x04
            0.0,  # par_x05
            0.5,  # par_x06
            5.0,  # par_x07
            0.0,  # par_x08
            0.0,  # par_x09
            0.0,  # par_x09+x10
            -5.0,  # 10**par_x11
            0.5,  # par_x12
            0.0,  # par_x13
            0.01,  # par_x13+x14
            0.005,  # par_x15
            -5.0,  # par_x16
            0.0,  # par_x17
            0.0,  # par_x18
            0.0,  # par_x19
            0.3,  # par_x20
            0.01,  # par_x21
            0.5,  # par_x22
            0.15,  # par_x23
            1.5,  # par_x24
            0.0,  # par_x24+x25
            -1.0,  # par_x26
            0.01,  # par_x27
            0.00001,  # par_x28
            0.0,  # par_x29
            0.0,  # par_x30
            -3.0,  # par_x31
            0.5,  # par_x32
            0.8,  # par_x33
            0.8,  # par_x34
            0.0,  # par_x35
            0.0,  # par_r01
            0.0,  # par_r02
            0.0,  # par_r03
            0.0,  # par_r04
            0.0,  # par_r05
            0.0,  # par_r06
            0.0,  # par_r07
            0.0,  # par_r08
        )
        high = (
            1.0,  # par_x01
            3.0,  # par_x02
            3.0,  # par_x03
            -1.0,  # 10**par_x04
            100.0,  # par_x05
            2.0,  # par_x06
            10.0,  # par_x07
            3.0,  # par_x08
            0.05,  # par_x09
            0.45,  # par_x09+x10
            -2.0,  # 10**par_x11
            2.0,  # par_x12
            0.1,  # par_x13
            0.3,  # par_x13+x14
            0.1,  # par_x15
            2.0,  # par_x16
            1.0,  # par_x17
            5.0,  # par_x18
            0.4,  # par_x19
            20.0,  # par_x20
            5.0,  # par_x21
            13.0,  # par_x22
            1.5,  # par_x23
            3.0,  # par_x24
            5.0,  # par_x24+x25
            1.0,  # par_x26
            0.2,  # par_x27
            0.02,  # par_x28
            0.5,  # par_x29
            2.0,  # par_x30
            3.0,  # par_x31
            4.0,  # par_x32
            1.2,  # par_x33
            1.2,  # par_x34
            0.02,  # par_x35
            1.0,  # par_r01
            1.0,  # par_r02
            1.0,  # par_r03
            1.0,  # par_r04
            1.0,  # par_r05
            1.0,  # par_r06
            1.0,  # par_r07
            1.0,  # par_r08
        )

        model.configure(
            get_local_testdata("ostrich-gr4j-cemaneige/OstRandomNumbers.txt")
        )

        model(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            hrus=(hru,),
            land_use_classes=(lu,),
            params=params,
            lowerBounds=low,
            upperBounds=high,
            algorithm="DDS",
            random_seed=0,
            max_iterations=10,
        )

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -1.47169, 4)
        # np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -1.51237, 4)   # true when x33 and x34 set to 1.0 and range [1.0,1.0]

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        expected_value = [
            2.400718e-02,
            1.516941e00,
            2.759658e00,
            -4.075778e00,
            3.483657e01,
            1.972741e00,
            7.476240e00,
            2.966089e00,
            2.035387e-03,
            3.379254e-01,
            -4.579974e00,
            6.867875e-01,
            8.920417e-02,
            1.681178e-01,
            8.827804e-02,
            -1.475022e00,
            4.722976e-01,
            4.528209e00,
            2.273521e-01,
            8.036873e00,
            3.461021e00,
            6.880423e00,
            1.312190e00,
            2.752630e00,
            1.515566e00,
            -1.499868e-01,
            6.437522e-02,
            1.013312e-02,
            1.089699e-01,
            1.462368e00,
            -1.620150e00,
            3.619720e00,
            1.130258e00,
            1.020023e00,
            1.622190e-02,
            7.319023e-02,
            1.081170e-01,
            1.222980e-01,
            4.622038e-01,
            4.863545e-02,
            3.800171e-01,
            9.620678e-01,
            7.091240e-01,
        ]
        np.testing.assert_almost_equal(
            opt_para,
            expected_value,
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            1.47169,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        blended = BLENDED()
        blended(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            hrus=(hru,),
            land_use_classes=(lu,),
            params=model.calibrated_params,
        )

        np.testing.assert_almost_equal(
            blended.diagnostics["DIAG_NASH_SUTCLIFFE"], d["DIAG_NASH_SUTCLIFFE"], 4
        )
