import numpy as np
import pandas as pd
import pytest

from ravenpy.utilities.regionalization import (
    distance,
    read_gauged_params,
    read_gauged_properties,
    regionalize,
)


class TestRegionalization:
    def test_full_example(self, symbolic_config):
        name, config = symbolic_config
        method = "SP_IDW"

        if name not in ["GR4JCN", "HMETS", "Mohyse"]:
            pytest.skip(f"Model {name} is not supported.")

        nash, params = read_gauged_params(name.upper())
        variables = ["latitude", "longitude", "area", "forest"]
        props = read_gauged_properties(variables)
        ungauged_props = {
            "latitude": 40.4848,
            "longitude": -103.3659,
            "area": 4250.6,
            "forest": 0.4,
        }

        qsim, ens = regionalize(
            config=config,
            method=method,
            nash=nash,
            params=params,
            props=props,
            target_props=ungauged_props,
            size=2,
            min_NSE=0.6,
        )

        assert qsim.max() > 1
        assert len(ens) == 2
        assert "members" in ens.dims
        assert "param" in ens.dims


class TestGreatCircle:
    def test_haversine_vector(self):
        places = dict()
        places["Ouagadougou"] = dict(latitude=12.366667, longitude=-1.533333)
        places["Nairobi"] = dict(latitude=-1.286389, longitude=36.817222)
        places["Cairo"] = dict(latitude=30.033333, longitude=31.233333)
        df = pd.DataFrame.from_dict(places).T
        algiers = pd.Series(dict(latitude=36.753889, longitude=3.058889))

        great_circle_vectors = distance(df, algiers)

        np.testing.assert_array_almost_equal(
            great_circle_vectors, [2750.249, 5478.388, 2709.166], 3
        )
