import datetime as dt
import tempfile

import numpy as np

from ravenpy.models import HRU, HYPR, HYPR_OST, LU
from ravenpy.utilities.testdata import get_local_testdata

from .common import _convert_2d

TS = get_local_testdata(
    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
)

hru = HYPR.HRU(
    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659, slope=0.01234
)

lu = LU("FOREST", impermeable_frac=0.0, forest_coverage=0.02345)


class TestHYPR:
    def test_simple(self):
        model = HYPR()
        params = (
            -1.856410e-01,  # par_x01
            2.92301100e00,  # par_x02
            3.1194200e-02,  # par_x03
            4.3982810e-01,  # par_x04
            4.6509760e-01,  # feed 10**par_x05; not par_x05=-3.324559E-01
            1.1770040e-01,  # feed 10**par_x06; not par_x06=-9.292220E-01
            1.31236800e01,  # par_x07
            4.0417950e-01,  # par_x08
            1.21225800e00,  # par_x09
            5.91273900e01,  # par_x10
            1.6612030e-01,  # par_x11
            4.10501500e00,  # par_x12
            8.2296110e-01,  # par_x13
            4.15635200e01,  # par_x14
            5.85111700e00,  # par_x15
            6.9090140e-01,  # par_x16
            9.2459950e-01,  # par_x17
            1.64358800e00,  # par_x18
            1.59920500e00,  # par_x19
            2.51938100e00,  # par_x20
            1.14820100e00,  # par_x21
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

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], 0.685188, 4)


class TestHYPR_OST:
    def test_simple(self):
        model = HYPR_OST()
        params = (
            -1.856410e-01,  # par_x01
            2.92301100e00,  # par_x02
            3.1194200e-02,  # par_x03
            4.3982810e-01,  # par_x04
            4.6509760e-01,  # feed 10**par_x05; not par_x05=-3.324559E-01
            1.1770040e-01,  # feed 10**par_x06; not par_x06=-9.292220E-01
            1.31236800e01,  # par_x07
            4.0417950e-01,  # par_x08
            1.21225800e00,  # par_x09
            5.91273900e01,  # par_x10
            1.6612030e-01,  # par_x11
            4.10501500e00,  # par_x12
            8.2296110e-01,  # par_x13
            4.15635200e01,  # par_x14
            5.85111700e00,  # par_x15
            6.9090140e-01,  # par_x16
            9.2459950e-01,  # par_x17
            1.64358800e00,  # par_x18
            1.59920500e00,  # par_x19
            2.51938100e00,  # par_x20
            1.14820100e00,  # par_x21
        )
        low = (
            -1.0,  # par_x01
            -3.0,  # par_x02
            0.00,  # par_x03
            0.30,  # par_x04
            -1.3,  # 10**par_x05
            -2.0,  # 10**par_x06
            0.00,  # par_x07
            0.10,  # par_x08
            0.40,  # par_x09
            0.00,  # par_x10
            0.00,  # par_x11
            0.00,  # par_x12
            0.00,  # par_x13
            0.00,  # par_x14
            0.01,  # par_x15
            0.00,  # par_x16
            0.00,  # par_x17
            1.50,  # par_x18
            0.00,  # par_x19
            0.00,  # par_x20
            0.80,  # par_x21
        )
        high = (
            1.0000,  # par_x01
            3.0000,  # par_x02
            0.8000,  # par_x03
            1.0000,  # par_x04
            0.3000,  # 10**par_x05
            0.0000,  # 10**par_x06
            30.000,  # par_x07
            0.8000,  # par_x08
            2.0000,  # par_x09
            100.00,  # par_x10
            0.5000,  # par_x11
            5.0000,  # par_x12
            1.0000,  # par_x13
            1000.0,  # par_x14
            6.0000,  # par_x15
            7.0000,  # par_x16
            8.0000,  # par_x17
            3.0000,  # par_x18
            5.0000,  # par_x19
            5.0000,  # par_x20
            1.2000,  # par_x21
        )

        model.configure(
            get_local_testdata("ostrich-gr4j-cemaneige/OstRandomNumbers.txt")
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
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

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -5.513390e-02, 4)

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        expected_value = [
            5.8879440e-01,  # par_x01
            -3.969576e-01,  # par_x02
            4.4449220e-01,  # par_x03
            4.5797310e-01,  # par_x04
            1.9331060e-01,  # 10**par_x05
            -1.4112350e00,  # 10**par_x06
            2.04891000e01,  # par_x07
            5.7430970e-01,  # par_x08
            6.7980080e-01,  # par_x09
            1.32665400e01,  # par_x10
            3.5348530e-01,  # par_x11
            2.7865950e-01,  # par_x12
            6.8593480e-01,  # par_x13
            5.05808700e02,  # par_x14
            2.4276360e-02,  # par_x15
            6.5476400e-01,  # par_x16
            7.0711200e-01,  # par_x17
            2.33075400e00,  # par_x18
            1.63560800e00,  # par_x19
            4.67153800e00,  # par_x20
            1.16302300e00,  # par_x21
        ]
        np.testing.assert_almost_equal(
            opt_para,
            expected_value,
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            5.513390e-02,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        hypr = HYPR()
        hypr(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            hrus=(hru,),
            land_use_classes=(lu,),
            params=model.calibrated_params,
        )

        np.testing.assert_almost_equal(
            hypr.diagnostics["DIAG_NASH_SUTCLIFFE"], d["DIAG_NASH_SUTCLIFFE"], 4
        )
