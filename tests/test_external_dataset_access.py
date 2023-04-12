import datetime as dt

import numpy as np
import pytest

from ravenpy.extractors.forecasts import get_CASPAR_dataset, get_ECCC_dataset


class TestGet:
    @pytest.mark.online
    def test_get_CASPAR_dataset(self):
        ds, _ = get_CASPAR_dataset("GEPS", dt.datetime(2018, 8, 31))

    @pytest.mark.online
    @pytest.mark.xfail(error=OSError, reason="Network may be unreliable")
    def test_get_ECCC_dataset(self):
        ds, _ = get_ECCC_dataset("GEPS")

        ns = np.datetime64("now") - ds.time.isel(time=0).values
        n_hours = ns / np.timedelta64(1, "h")

        assert n_hours <= 36
