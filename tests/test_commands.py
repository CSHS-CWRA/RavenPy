from ravenpy.config import commands as c
from ravenpy.utilities.testdata import get_local_testdata


class TestDataCommand:
    def test_infer_scale_and_offset(self):
        ts = get_local_testdata("era5/tas_pr_20180101-20180108.nc")
        pr = c.DataCommand(data_type="PRECIP", file_name_nc=ts, var_name_nc="pr")

        assert pr.scale is None

        pr.infer_scale_and_offset()
        assert pr.scale == 24000.0
        assert pr.offset == 0

        tas = c.DataCommand(
            data_type="TEMP_DAILY_AVE", file_name_nc=ts, var_name_nc="tas"
        )
        tas.infer_scale_and_offset()
        assert tas.scale == 1
        assert tas.offset == -273.15
