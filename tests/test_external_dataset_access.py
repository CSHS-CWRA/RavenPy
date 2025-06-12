import datetime as dt
import urllib.error
from pathlib import Path

import numpy as np
import pytest
import xarray

from ravenpy.extractors.forecasts import get_CASPAR_dataset, get_ECCC_dataset
from ravenpy.testing.utils import default_testdata_cache, open_dataset, yangtze


@pytest.mark.online
class TestGet:
    def test_get_CASPAR_dataset(self):
        ds, _ = get_CASPAR_dataset("GEPS", dt.datetime(2018, 8, 31))

    @pytest.mark.xfail(error=OSError, reason="Network may be unreliable", strict=False)
    def test_get_ECCC_dataset(self):
        ds, _ = get_ECCC_dataset("GEPS")

        ns = np.datetime64("now") - ds.time.isel(time=0).values
        n_hours = ns / np.timedelta64(1, "h")

        assert n_hours <= 36


@pytest.mark.online
class TestRemoteFileAccess:
    dap_url = "http://test.opendap.org:80/opendap/data/nc/"
    git_url = "https://github.com/Ouranosinc/raven-testdata"
    branch = "main"

    @pytest.mark.xfail(
        raises=urllib.error.URLError,
        reason="Get file is API rate limited",
        strict=False,
    )
    def test_get_file_default_cache(self):
        file = yangtze(branch=self.branch).fetch(
            fname="ostrich-hbvec/raven-hbvec-salmon.rvi"
        )

        assert Path(default_testdata_cache).exists()
        assert Path(file).is_file()
        with Path(file).open() as f:
            header = f.read()
            assert ":FileType          rvi ASCII Raven 2.8.2" in header

    @pytest.mark.xfail(
        raises=urllib.error.HTTPError,
        reason="Transitory failures expected",
        strict=False,
    )
    def test_open_dataset(self):
        ds = open_dataset(
            name="raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
            yangtze_kwargs={"branch": self.branch, "force_download": True},
        )

        assert (
            Path(default_testdata_cache)
            .joinpath(
                self.branch,
                "raven-gr4j-cemaneige",
                "Salmon-River-Near-Prince-George_meteo_daily.nc",
            )
            .exists()
        )
        assert isinstance(ds, xarray.Dataset)
