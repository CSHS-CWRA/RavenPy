import xarray as xr

from ravenpy.utilities.coords import infer_scale_and_offset
from ravenpy.utilities.testdata import get_local_testdata


def test_infer_scale_and_offset():
    # ERA5 precip and tas
    ts = get_local_testdata("era5/tas_pr_20180101-20180108.nc")
    with xr.open_dataset(ts) as ds:
        p = infer_scale_and_offset(ds["pr"], "PRECIP")
        assert p == (24000, 0)

        p = infer_scale_and_offset(ds["tas"], "TEMP_DAILY_AVE")
        assert p == (1, -273.15)

    # CMIP precip in kg / m^2 / s
    ts2 = get_local_testdata(
        "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp45_nex-gddp_1971-1972_subset.nc"
    )
    with xr.open_dataset(ts2) as ds:
        p = infer_scale_and_offset(ds["pr"], "PRECIP")
        assert p == (86400, 0)

    # VIC streamflow inputs in mm accumulated over 6 hours
    fn = get_local_testdata("raven-routing-sample/VIC_streaminputs.nc")
    with xr.open_dataset(fn) as ds:
        p = infer_scale_and_offset(ds.Streaminputs, "PRECIP")
        assert p == (4, 0)
