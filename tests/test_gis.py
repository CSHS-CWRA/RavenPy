import pytest
from ravenpy.utilities import gis


@pytest.mark.skip(
    reason="The hydrobasins_domains dataset associated to these tests will be made eventually available through a `raven-testdata` repository."
)
class TestSelect_hybas_domain:
    def test_na(self):
        bbox = (-68.0, 50.0, -68.0, 50.0)
        dom = gis.select_hybas_domain(bbox)
        assert dom == "na"

    def test_ar(self):
        bbox = (-114.65, 61.35, -114.65, 61.35)
        dom = gis.select_hybas_domain(bbox)
        assert dom == "ar"
