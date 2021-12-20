import datetime as dt

import numpy as np

from ravenpy.config.commands import LU
from ravenpy.models import SACSMA, SACSMA_OST
from ravenpy.utilities.testdata import get_local_testdata

TS = get_local_testdata(
    "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
)

hru = SACSMA.HRU(
    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659, slope=0.01234
)

lu = LU("FOREST", impermeable_frac=0.0, forest_coverage=0.02345)


class TestSACSMA:
    def test_simple(self):
        model = SACSMA()
        params = SACSMA.Params(
            0.0100000,  # feed 10**par_x01; ; not par_x1=???
            0.0500000,  # feed 10**par_x02; ; not par_x2=???
            0.3000000,  # feed 10**par_x03; ; not par_x3=???
            0.0500000,  # par_x04
            0.0500000,  # par_x05
            0.1300000,  # par_x06
            0.0250000,  # par_x07
            0.0600000,  # par_x08
            0.0600000,  # par_x09
            1.0000000,  # par_x10
            40.000000,  # par_x11
            0.0000000,  # feed par_x12/(1+par_x12) par_x12; not par_x12=???
            0.0000000,  # par_x13
            0.1000000,  # par_x14
            0.0000000,  # par_x15
            0.0100000,  # par_x16
            1.5000000,  # par_x17
            0.4827523,  # par_x18
            4.0998200,  # par_x19
            1.0000000,  # par_x20
            1.0000000,  # par_x21
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

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], -0.0382907, 6)


class TestSACSMA_OST:
    def test_simple(self):
        model = SACSMA_OST()
        params = SACSMA.Params(
            0.0100000,  # feed 10**par_x01; ; not par_x1=???
            0.0500000,  # feed 10**par_x02; ; not par_x2=???
            0.3000000,  # feed 10**par_x03; ; not par_x3=???
            0.0500000,  # par_x04
            0.0500000,  # par_x05
            0.1300000,  # par_x06
            0.0250000,  # par_x07
            0.0600000,  # par_x08
            0.0600000,  # par_x09
            1.0000000,  # par_x10
            40.000000,  # par_x11
            0.0000000,  # feed par_x12/(1+par_x12) par_x12; not par_x12=???
            0.0000000,  # par_x13
            0.1000000,  # par_x14
            0.0000000,  # par_x15
            0.0100000,  # par_x16
            1.5000000,  # par_x17
            0.4827523,  # par_x18
            4.0998200,  # par_x19
            1.0000000,  # par_x20
            1.0000000,  # par_x21
        )
        low = SACSMA.Params(
            -3.00000000,  # 10**val = 0.001
            -1.52287874,  # 10**val = 0.03
            -0.69897000,  # 10**val = 0.2
            0.025000000,  # par_x04
            0.010000000,  # par_x05
            0.075000000,  # par_x06
            0.015000000,  # par_x07
            0.040000000,  # par_x08
            0.000000000,  # par_x09
            0.000000000,  # par_x10
            0.000000000,  # par_x11
            0.000000000,  # par_x12/(1+par_x12)
            0.000000000,  # par_x13
            0.000000000,  # par_x14
            0.000000000,  # par_x15
            0.000000000,  # par_x16
            0.000000000,  # par_x17
            0.300000000,  # par_x18
            0.010000000,  # par_x19
            0.800000000,  # par_x20
            0.800000000,  # par_x21
        )
        high = SACSMA.Params(
            -1.82390874,  # 10**val = 0.015
            -0.69897000,  # 10**val = 0.2
            -0.30102999,  # 10**val = 0.5
            0.125000000,  # par_x04
            0.075000000,  # par_x05
            0.300000000,  # par_x06
            0.300000000,  # par_x07
            0.600000000,  # par_x08
            0.500000000,  # par_x09
            3.000000000,  # par_x10
            80.00000000,  # par_x11
            0.800000000,  # par_x12/(1+par_x12)
            0.050000000,  # par_x13
            0.200000000,  # par_x14
            0.100000000,  # par_x15
            0.400000000,  # par_x16
            8.000000000,  # par_x17
            20.00000000,  # par_x18
            5.000000000,  # par_x19
            1.200000000,  # par_x20
            1.200000000,  # par_x21
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

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], 0.532016, 4)

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        expected_value = [
            -2.065716e00,
            -1.165434e00,
            -4.778685e-01,
            4.756758e-02,
            7.066574e-02,
            1.412360e-01,
            2.096465e-01,
            4.194477e-01,
            7.903211e-02,
            3.979961e-01,
            5.655765e01,
            4.458551e-02,
            3.429674e-02,
            1.011617e-01,
            2.383366e-04,
            3.741509e-02,
            7.071120e-01,
            1.121058e01,
            1.642337e00,
            1.173723e00,
            1.163023e00,
        ]
        np.testing.assert_almost_equal(
            opt_para,
            expected_value,
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            -0.532016,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        sacsma = SACSMA()
        sacsma(
            TS,
            start_date=dt.datetime(1954, 1, 1),
            duration=208,
            hrus=(hru,),
            land_use_classes=(lu,),
            params=model.calibrated_params,
        )

        np.testing.assert_almost_equal(
            sacsma.diagnostics["DIAG_NASH_SUTCLIFFE"], d["DIAG_NASH_SUTCLIFFE"], 4
        )
