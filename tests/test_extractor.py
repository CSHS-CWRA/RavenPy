import pytest

from ravenpy.config.emulators import BasicRoute


@pytest.mark.gis
class TestExtractor:
    routing_product = pytest.importorskip("ravenpy.extractors.routing_product")

    def test_basinmaker_extractor(self, tmp_path, yangtze):

        routing_product_shp_path = yangtze.fetch("basinmaker/drainage_region_0175_v2-1/finalcat_info_v2-1.zip")
        df = self.routing_product.open_shapefile(
            routing_product_shp_path,
        )
        rvh_extractor = self.routing_product.BasinMakerExtractor(
            df=df,
            hru_aspect_convention="ArcGIS",
            routing_product_version="2.1",
        )
        rvh_config = rvh_extractor.extract(hru_from_sb=True)

        # Create lists of values to check
        bedslope_list = [item["bed_slope"] for item in rvh_config["channel_profile"]]
        mannings_list = [value for d in rvh_config["channel_profile"] for value in [t[1] for t in d["roughness_zones"]]]
        reach_length_list = [item["reach_length"] for item in rvh_config["sub_basins"]]

        rvh_config.pop("channel_profile")

        config = BasicRoute(**rvh_config)
        config.write_rv(tmp_path, modelname="routing")

        # Checks that the bedslope, Manning and reach length values are non negative numbers
        assert all(isinstance(x, (int, float)) for x in bedslope_list) is True
        assert any(x < 0 for x in bedslope_list) is False
        assert all(isinstance(y, (int, float)) for y in mannings_list) is True
        assert any(y < 0 for y in mannings_list) is False
        assert all(isinstance(z, (int, float)) for z in bedslope_list) is True
        assert any(z < 0 for z in reach_length_list) is False
