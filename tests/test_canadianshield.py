import datetime as dt
import tempfile

import numpy as np
import pytest

from ravenpy.config import ConfigError
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
        model = CANADIANSHIELD()
        model.config.rvh.hrus = (
            CANADIANSHIELD.HRU_ORGANIC(**hru_default_values),
            CANADIANSHIELD.HRU_BEDROCK(**hru_default_values),
        )
        params = (
            4.72304300e-01,  # par_x01
            8.16392200e-01,  # par_x02
            9.86197600e-02,  # par_x03
            3.92699900e-03,  # par_x04
            4.69073600e-02,  # par_x05
            4.95528400e-01,  # feed par_x05+x06; not par_x06=4.4862100E-01
            6.803492000e00,  # par_x07
            4.33050200e-03,  # feed 10**par_x08; not par_x8=-2.363462E+00
            1.01425900e-05,  # feed 10**par_x09; not par_x9=-4.993851E+00
            1.823470000e00,  # par_x10
            5.12215400e-01,  # par_x11
            9.017555000e00,  # par_x12
            3.077103000e01,  # par_x13
            5.094095000e01,  # par_x14
            1.69422700e-01,  # par_x15
            8.23412200e-02,  # par_x16
            2.34595300e-01,  # feed par_x16+x17; not par_x17=1.5225400E-01
            7.30904000e-02,  # par_x18
            1.284052000e00,  # par_x19
            3.653415000e00,  # par_x20
            2.306515000e01,  # par_x21
            2.402183000e00,  # par_x22
            2.522095000e00,  # par_x23
            5.80344900e-01,  # par_x24
            1.614157000e00,  # par_x25
            6.031781000e00,  # par_x26
            3.11129800e-01,  # par_x27
            6.71695100e-02,  # par_x28
            5.83759500e-05,  # par_x29
            9.824723000e00,  # par_x30
            9.00747600e-01,  # par_x31
            8.04057300e-01,  # par_x32
            1.179003000e00,  # par_x33
            7.98001300e-01,  # par_x34
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            land_use_classes=(lu,),
            params=params,
            suppress_output=True,
        )

        d = model.diagnostics

        np.testing.assert_almost_equal(d["DIAG_NASH_SUTCLIFFE"], 0.39602, 4)

    def test_bad_config(self):
        model = CANADIANSHIELD()
        params = (
            4.72304300e-01,  # par_x01
            8.16392200e-01,  # par_x02
            9.86197600e-02,  # par_x03
            3.92699900e-03,  # par_x04
            4.69073600e-02,  # par_x05
            4.95528400e-01,  # feed par_x05+x06; not par_x06=4.4862100E-01
            6.803492000e00,  # par_x07
            4.33050200e-03,  # feed 10**par_x08; not par_x8=-2.363462E+00
            1.01425900e-05,  # feed 10**par_x09; not par_x9=-4.993851E+00
            1.823470000e00,  # par_x10
            5.12215400e-01,  # par_x11
            9.017555000e00,  # par_x12
            3.077103000e01,  # par_x13
            5.094095000e01,  # par_x14
            1.69422700e-01,  # par_x15
            8.23412200e-02,  # par_x16
            2.34595300e-01,  # feed par_x16+x17; not par_x17=1.5225400E-01
            7.30904000e-02,  # par_x18
            1.284052000e00,  # par_x19
            3.653415000e00,  # par_x20
            2.306515000e01,  # par_x21
            2.402183000e00,  # par_x22
            2.522095000e00,  # par_x23
            5.80344900e-01,  # par_x24
            1.614157000e00,  # par_x25
            6.031781000e00,  # par_x26
            3.11129800e-01,  # par_x27
            6.71695100e-02,  # par_x28
            5.83759500e-05,  # par_x29
            9.824723000e00,  # par_x30
            9.00747600e-01,  # par_x31
            8.04057300e-01,  # par_x32
            1.179003000e00,  # par_x33
            7.98001300e-01,  # par_x34
        )

        # Must be 2 HRUs
        with pytest.raises(ConfigError) as _:
            model.config.rvh.hrus = (CANADIANSHIELD.HRU_ORGANIC(**hru_default_values),)
            model(
                TS,
                start_date=dt.datetime(2000, 1, 1),
                end_date=dt.datetime(2002, 1, 1),
                land_use_classes=(lu,),
                params=params,
                suppress_output=True,
            )

        # The 2 HRUs must have the same area
        with pytest.raises(ConfigError) as _:
            model.config.rvh.hrus = (
                CANADIANSHIELD.HRU_ORGANIC(area=100),
                CANADIANSHIELD.HRU_BEDROCK(area=200),
            )
            model(
                TS,
                start_date=dt.datetime(2000, 1, 1),
                end_date=dt.datetime(2002, 1, 1),
                land_use_classes=(lu,),
                params=params,
                suppress_output=True,
            )


class TestCANADIANSHIELD_OST:
    def test_simple(self):
        model = CANADIANSHIELD_OST()
        model.config.rvh.hrus = (
            CANADIANSHIELD.HRU_ORGANIC(**hru_default_values),
            CANADIANSHIELD.HRU_BEDROCK(**hru_default_values),
        )
        params = (
            4.72304300e-01,  # par_x01
            8.16392200e-01,  # par_x02
            9.86197600e-02,  # par_x03
            3.92699900e-03,  # par_x04
            4.69073600e-02,  # par_x05
            4.95528400e-01,  # feed par_x05+x06; not par_x06=4.4862100E-01
            6.803492000e00,  # par_x07
            4.33050200e-03,  # feed 10**par_x08; not par_x8=-2.363462E+00
            1.01425900e-05,  # feed 10**par_x09; not par_x9=-4.993851E+00
            1.823470000e00,  # par_x10
            5.12215400e-01,  # par_x11
            9.017555000e00,  # par_x12
            3.077103000e01,  # par_x13
            5.094095000e01,  # par_x14
            1.69422700e-01,  # par_x15
            8.23412200e-02,  # par_x16
            2.34595300e-01,  # feed par_x16+x17; not par_x17=1.5225400E-01
            7.30904000e-02,  # par_x18
            1.284052000e00,  # par_x19
            3.653415000e00,  # par_x20
            2.306515000e01,  # par_x21
            2.402183000e00,  # par_x22
            2.522095000e00,  # par_x23
            5.80344900e-01,  # par_x24
            1.614157000e00,  # par_x25
            6.031781000e00,  # par_x26
            3.11129800e-01,  # par_x27
            6.71695100e-02,  # par_x28
            5.83759500e-05,  # par_x29
            9.824723000e00,  # par_x30
            9.00747600e-01,  # par_x31
            8.04057300e-01,  # par_x32
            1.179003000e00,  # par_x33
            7.98001300e-01,  # par_x34
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

        model.configure(
            get_local_testdata("ostrich-gr4j-cemaneige/OstRandomNumbers.txt")
        )

        model(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
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

        opt_para = model.optimized_parameters
        opt_func = model.obj_func

        expected_value = [
            3.1462940e-01,
            1.46843200e00,
            3.1906490e-01,
            1.25061900e00,
            3.6118660e-02,
            4.2020140e-01,
            6.55302700e00,
            -1.7296660e00,
            -4.5009800e00,
            1.63844700e00,
            1.87916000e00,
            3.54916500e00,
            4.88600200e01,
            9.03863300e01,
            9.2422170e-02,
            2.0133650e-02,
            2.9473000e-01,
            5.9477830e-02,
            1.27923400e00,
            5.6596230e-01,
            3.76718000e02,
            7.0004300e-01,
            1.46996300e00,
            7.4085920e-01,
            5.75679100e00,
            1.49193800e01,
            7.5535240e-01,
            9.4459520e-02,
            1.8634570e-01,
            6.6940050e-01,
            3.92734700e00,
            1.00278400e00,
            1.00417400e00,
            7.8267860e-01,
        ]
        np.testing.assert_almost_equal(
            opt_para,
            expected_value,
            4,
            err_msg="calibrated parameter set is not matching expected value",
        )
        np.testing.assert_almost_equal(
            opt_func,
            4.790550e-01,
            4,
            err_msg="calibrated NSE is not matching expected value",
        )

        canadianshield = CANADIANSHIELD()
        canadianshield.config.rvh.hrus = (
            CANADIANSHIELD.HRU_ORGANIC(**hru_default_values),
            CANADIANSHIELD.HRU_BEDROCK(**hru_default_values),
        )
        canadianshield(
            TS,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            land_use_classes=(lu,),
            params=model.calibrated_params,
            suppress_output=True,
        )

        np.testing.assert_almost_equal(
            canadianshield.diagnostics["DIAG_NASH_SUTCLIFFE"],
            d["DIAG_NASH_SUTCLIFFE"],
            4,
        )
