import xarray as xr
import xclim.sdba
import xclim.sdba as sdba
from xclim.core.calendar import convert_calendar

from ravenpy.utilities.testdata import get_local_testdata


class TestBiasCorrect:
    def test_bias_correction(self):

        ds_fut_sub = xr.open_dataset(
            get_local_testdata(
                "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp85_nex-gddp_2070-2071_subset.nc",
            )
        )
        ds_ref_sub = xr.open_dataset(
            get_local_testdata(
                "cmip5/nasa_nex-gddp-1.0_day_inmcm4_historical+rcp45_nex-gddp_1971-1972_subset.nc",
            )
        )
        ds_ref_sub = convert_calendar(ds_ref_sub, "noleap")

        ds_his_sub = xr.open_dataset(
            get_local_testdata("nrcan/NRCAN_1971-1972_subset.nc")
        )
        ds_his_sub = convert_calendar(ds_his_sub, "noleap")
        group = xclim.sdba.Grouper("time.month")
        # Train the model to find the correction factors
        Adj = sdba.DetrendedQuantileMapping.train(
            ref=ds_ref_sub["pr"],
            hist=ds_his_sub["pr"],
            nquantiles=50,
            kind="+",
            group=group,
        )

        # Apply the factors to the future data to bias-correct
        Adj.adjust(ds_fut_sub["pr"], interp="linear")

        # Repeat for temperature max
        Adj = sdba.DetrendedQuantileMapping.train(
            ref=ds_ref_sub["tasmax"],
            hist=ds_his_sub["tasmax"],
            nquantiles=50,
            kind="+",
            group=group,
        )

        # Apply the factors to the future data to bias-correct
        Adj.adjust(ds_fut_sub["tasmax"], interp="linear")

        # Repeat for tasmin
        Adj = sdba.DetrendedQuantileMapping.train(
            ref=ds_ref_sub["tasmin"],
            hist=ds_his_sub["tasmin"],
            nquantiles=50,
            kind="+",
            group=group,
        )

        Adj.adjust(ds_fut_sub["tasmin"], interp="linear")

        # TODO: Add numerical check
