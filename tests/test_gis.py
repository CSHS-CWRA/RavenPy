import tempfile

import fiona
import numpy as np
import pytest
from shapely.geometry import shape

from ravenpy.utilities import gis


# @pytest.mark.skip(
#     reason="The hydrobasins_domains dataset associated to these tests will be made eventually available through a "
#            "`raven-testdata` repository."
# )
class TestHydroBASINS:
    def test_select_hybas_na_domain(self):
        bbox = (-68.0, 50.0, -68.0, 50.0)
        dom = gis.select_hybas_domain(bbox)
        assert dom == "na"

    def test_select_hybas_ar_domain(self):
        bbox = (-114.65, 61.35, -114.65, 61.35)
        dom = gis.select_hybas_domain(bbox)
        assert dom == "ar"

    def test_hydrobasins_aggregate(self):
        pass

    def test_hydrobasins_upstream_ids(self):
        pass

    def test_get_hydrobasins_location_wfs(self):
        lake_winnipeg = (-98.03575958286369, 52.88238524279493)
        feature = gis.get_hydrobasins_location_wfs(
            coordinates=lake_winnipeg * 2, lakes=True, domain="na"
        )

        gml = tempfile.NamedTemporaryFile(suffix=".gml", delete=False)
        with open(gml.name, "wb") as f:
            f.write(feature)

        with fiona.open(gml.name) as src:
            iterable = iter(src)
            feat = next(iterable)
            geom = shape(feat["geometry"])
            assert geom.bounds == (-99.2731, 50.3603, -96.2578, 53.8705)
            np.testing.assert_almost_equal(geom.area, 3.2530867)


class TestHydroRouting:
    @pytest.mark.skip(reason="This function doesn't exist yet.")
    def test_hydro_routing_aggregate(self):
        pass

    def test_hydro_routing_upstream_ids(self):
        pass


class TestWFS:
    def test_get_location_wfs(self):
        pass

    def test_get_feature_attributes_wfs(self):
        pass


class TestWCS:
    def test_get_raster_wcs(self):
        pass


class TestGIS:
    def test_get_bbox(self):
        pass

    def test_feature_contains(self):
        pass
