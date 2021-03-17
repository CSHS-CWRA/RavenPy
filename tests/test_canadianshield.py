import datetime as dt
import tempfile

import numpy as np

from ravenpy.models import CANADIANSHIELD, CANADIANSHIELD_OST, HRU, LU
from ravenpy.utilities.testdata import get_local_testdata

from .common import _convert_2d

TS = get_local_testdata(
    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
)

lu = LU("FOREST", impermeable_frac=0.0, forest_coverage=0.02345)


hru_default_values = dict(
    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659, slope=0.01234
)


class TestCANADIANSHIELD:
    def test_simple(self):
        model = CANADIANSHIELD()  # this is running __init__ in emulator
        model.rvh.hrus = (
            CANADIANSHIELD.HRU_ORGANIC(**hru_default_values),
            CANADIANSHIELD.HRU_BEDROCK(**hru_default_values),
        )
        params = (
            4.7230430e-01,  # par_x01
            8.1639220e-01,  # par_x02
            9.8619760e-02,  # par_x03
            3.9269990e-03,  # par_x04
            4.6907360e-02,  # par_x05
            4.9552840e-01,  # feed par_x05+x06; not par_x06=4.4862100E-01  4.955284E-01
            6.8034920e00,  # par_x07
            4.3305020e-03,  # feed 10**par_x08; not par_x8=-2.363462E+00   4.330502E-03
            1.0142590e-05,  # feed 10**par_x09; not par_x9=-4.993851E+00   1.014259E-05
            1.8234700e00,  # par_x10
            5.1221540e-01,  # par_x11
            9.0175550e00,  # par_x12
            3.0771030e01,  # par_x13
            5.0940950e01,  # par_x14
            1.6942270e-01,  # par_x15
            8.2341220e-02,  # par_x16
            2.3459530e-01,  # feed par_x16+x17; not par_x17=1.5225400E-01  2.345953E-01
            7.3090400e-02,  # par_x18
            1.2840520e00,  # par_x19
            3.6534150e00,  # par_x20
            2.3065150e01,  # par_x21
            2.4021830e00,  # par_x22
            2.5220950e00,  # par_x23
            5.8034490e-01,  # par_x24
            1.6141570e00,  # par_x25
            6.0317810e00,  # par_x26
            3.1112980e-01,  # par_x27
            6.7169510e-02,  # par_x28
            5.8375950e-05,  # par_x29
            9.8247230e00,  # par_x30
            9.0074760e-01,  # par_x31
            8.0405730e-01,  # par_x32
            1.1790030e00,  # par_x33
            7.9800130e-01,  # par_x34
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            # hrus=(hru,),
            land_use_classes=(lu,),
            params=params,
            suppress_output=True,
        )

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], 0.39602, 4)


class TestCANADIANSHIELD_OST:
    def test_simple(self):
        model = CANADIANSHIELD_OST("/tmp/pavics-julie/")
        model.rvh.hrus = (
            CANADIANSHIELD.HRU_ORGANIC(**hru_default_values),
            CANADIANSHIELD.HRU_BEDROCK(**hru_default_values),
        )
        params = (
            4.7230430e-01,  # par_x01
            8.1639220e-01,  # par_x02
            9.8619760e-02,  # par_x03
            3.9269990e-03,  # par_x04
            4.6907360e-02,  # par_x05
            4.9552840e-01,  # feed par_x05+x06; not par_x06=4.4862100E-01
            6.8034920e00,  # par_x07
            4.3305020e-03,  # feed 10**par_x08; not par_x8=-2.363462E+00
            1.0142590e-05,  # feed 10**par_x09; not par_x9=-4.993851E+00
            1.8234700e00,  # par_x10
            5.1221540e-01,  # par_x11
            9.0175550e00,  # par_x12
            3.0771030e01,  # par_x13
            5.0940950e01,  # par_x14
            1.6942270e-01,  # par_x15
            8.2341220e-02,  # par_x16
            2.3459530e-01,  # feed par_x16+x17; not par_x17=1.5225400E-01
            7.3090400e-02,  # par_x18
            1.2840520e00,  # par_x19
            3.6534150e00,  # par_x20
            2.3065150e01,  # par_x21
            2.4021830e00,  # par_x22
            2.5220950e00,  # par_x23
            5.8034490e-01,  # par_x24
            1.6141570e00,  # par_x25
            6.0317810e00,  # par_x26
            3.1112980e-01,  # par_x27
            6.7169510e-02,  # par_x28
            5.8375950e-05,  # par_x29
            9.8247230e00,  # par_x30
            9.0074760e-01,  # par_x31
            8.0405730e-01,  # par_x32
            1.1790030e00,  # par_x33
            7.9800130e-01,  # par_x34
        )
        low = (
            0.010,  # par_x01
            0.010,  # par_x02
            0.010,  # par_x03
            0.000,  # par_x04
            0.000,  # par_x05
            0.050,  # feed par_x05+x06; not par_x06=4.4862100E-01
            0.000,  # par_x07
            -5.00,  # feed 10**par_x08; not par_x8=-2.363462E+00
            -5.00,  # feed 10**par_x09; not par_x9=-4.993851E+00
            0.500,  # par_x10
            0.500,  # par_x11
            0.000,  # par_x12
            0.000,  # par_x13
            0.000,  # par_x14
            0.000,  # par_x15
            0.000,  # par_x16
            0.010,  # feed par_x16+x17; not par_x17=1.5225400E-01
            0.005,  # par_x18
            -3.00,  # par_x19
            0.500,  # par_x20
            5.000,  # par_x21
            0.000,  # par_x22
            0.000,  # par_x23
            -1.00,  # par_x24
            0.000,  # par_x25
            0.000,  # par_x26
            0.000,  # par_x27
            0.000,  # par_x28
            0.000,  # par_x29
            0.000,  # par_x30
            0.000,  # par_x31
            0.800,  # par_x32
            0.800,  # par_x33
            0.000,  # par_x34
        )
        high = (
            0.500,  # par_x01
            2.000,  # par_x02
            3.000,  # par_x03
            3.000,  # par_x04
            0.050,  # par_x05
            0.450,  # feed par_x05+x06; not par_x06=4.4862100E-01  4.955284E-01
            7.000,  # par_x07
            -1.00,  # feed 10**par_x08; not par_x8=-2.363462E+00   4.330502E-03
            -1.00,  # feed 10**par_x09; not par_x9=-4.993851E+00   1.014259E-05
            2.000,  # par_x10
            2.000,  # par_x11
            100.0,  # par_x12
            100.0,  # par_x13
            100.0,  # par_x14
            0.400,  # par_x15
            0.100,  # par_x16
            0.300,  # feed par_x16+x17; not par_x17=1.5225400E-01  2.345953E-01
            0.100,  # par_x18
            3.000,  # par_x19
            4.000,  # par_x20
            500.0,  # par_x21
            5.000,  # par_x22
            5.000,  # par_x23
            1.000,  # par_x24
            8.000,  # par_x25
            20.00,  # par_x26
            1.500,  # par_x27
            0.200,  # par_x28
            0.200,  # par_x29
            10.00,  # par_x30
            10.00,  # par_x31
            1.200,  # par_x32
            1.200,  # par_x33
            1.000,  # par_x34
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            # hrus=(hru,),
            land_use_classes=(lu,),
            params=params,
            lowerBounds=low,
            upperBounds=high,
            algorithm="DDS",
            random_seed=0,
            max_iterations=10,
        )

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -4.790550e-01, 4)


#         opt_para = model.optimized_parameters
#         opt_func = model.obj_func

#         expected_value = [
#             2.400718E-02,
#             1.516941E+00,
#             2.759658E+00,
#             -4.075778E+00,
#             3.483657E+01,
#             1.972741E+00,
#             7.476240E+00,
#             2.966089E+00,
#             2.035387E-03,
#             3.379254E-01,
#             -4.579974E+00,
#             6.867875E-01,
#             8.920417E-02,
#             1.681178E-01,
#             8.827804E-02,
#             -1.475022E+00,
#             4.722976E-01,
#             4.528209E+00,
#             2.273521E-01,
#             8.036873E+00,
#             3.461021E+00,
#             6.880423E+00,
#             1.312190E+00,
#             2.752630E+00,
#             1.515566E+00,
#             -1.499868E-01,
#             6.437522E-02,
#             1.013312E-02,
#             1.089699E-01,
#             1.462368E+00,
#             -1.620150E+00,
#             3.619720E+00,
#             1.130258E+00,
#             1.020023E+00,
#             1.622190E-02,
#             7.319023E-02,
#             1.081170E-01,
#             1.222980E-01,
#             4.622038E-01,
#             4.863545E-02,
#             3.800171E-01,
#             9.620678E-01,
#             7.091240E-01
#         ]
#         np.testing.assert_almost_equal(
#             opt_para,
#             expected_value,
#             4,
#             err_msg="calibrated parameter set is not matching expected value",
#         )
#         np.testing.assert_almost_equal(
#             opt_func,
#             1.47169,
#             4,
#             err_msg="calibrated NSE is not matching expected value",
#         )

#         canadianshield = CANADIANSHIELD()
#         canadianshield(
#             TS,
#             start_date=dt.datetime(1954, 1, 1),
#             duration=208,
#             hrus=(hru,),
#             land_use_classes=(lu,),
#             params=model.calibrated_params,
#         )

#         np.testing.assert_almost_equal(
#             canadianshield.diagnostics["DIAG_NASH_SUTCLIFFE"], d["DIAG_NASH_SUTCLIFFE"], 4
#         )
