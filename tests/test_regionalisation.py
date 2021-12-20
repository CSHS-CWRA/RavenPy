import datetime as dt

import numpy as np
import pandas as pd

from ravenpy.models import GR4JCN
from ravenpy.utilities import regionalization as reg
from ravenpy.utilities.testdata import get_local_testdata


class TestRegionalization:
    def test_full_example(self):
        ts = get_local_testdata(
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
        )
        model = "GR4JCN"
        nash, params = reg.read_gauged_params(model)
        variables = ["latitude", "longitude", "area", "forest"]
        props = reg.read_gauged_properties(variables)
        ungauged_props = {
            "latitude": 40.4848,
            "longitude": -103.3659,
            "area": 4250.6,
            "forest": 0.4,
        }

        hrus = (
            GR4JCN.LandHRU(
                area=4250.6, elevation=843.0, latitude=40.4848, longitude=-103.3659
            ),
        )

        qsim, ens = reg.regionalize(
            "SP_IDW",
            model,
            nash,
            params,
            props,
            ungauged_props,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 1, 1),
            hrus=hrus,
            longitude=-103.3659,
            min_NSE=0.6,
            size=2,
            ts=ts,
        )

        assert qsim.max() > 1
        assert len(ens) == 2
        assert "realization" in ens.dims
        assert "param" in ens.dims


class TestGreatCircle:
    def test_haversine_vector(self):
        places = dict()
        places["Ouagadougou"] = dict(latitude=12.366667, longitude=-1.533333)
        places["Nairobi"] = dict(latitude=-1.286389, longitude=36.817222)
        places["Cairo"] = dict(latitude=30.033333, longitude=31.233333)
        df = pd.DataFrame.from_dict(places).T
        algiers = pd.Series(dict(latitude=36.753889, longitude=3.058889))

        great_circle_vectors = reg.distance(df, algiers)

        np.testing.assert_array_almost_equal(
            great_circle_vectors, [2750.249, 5478.388, 2709.166], 3
        )
