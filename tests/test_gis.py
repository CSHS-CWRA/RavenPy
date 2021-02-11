import pytest
from ravenpy.utilities import gis


@pytest.mark.skip(
    reason="The hydrobasins_domains dataset associated to these tests will be made eventually available through a "
           "`raven-testdata` repository."
)
class TestSelectHydroBASINSDomain:
    def test_na(self):
        bbox = (-68.0, 50.0, -68.0, 50.0)
        dom = gis.select_hybas_domain(bbox)
        assert dom == "na"

    def test_ar(self):
        bbox = (-114.65, 61.35, -114.65, 61.35)
        dom = gis.select_hybas_domain(bbox)
        assert dom == "ar"


class TestHydroBASINSAggregate:
    pass


class TestHydroBASINSUpstreamIDs:
    pass


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
