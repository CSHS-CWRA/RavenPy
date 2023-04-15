from ravenpy.config.emulators import BasicRoute
from ravenpy.extractors.routing_product import BasinMakerExtractor, open_shapefile


def test_basinmaker_extractor(get_local_testdata, tmp_path):
    routing_product_shp_path = get_local_testdata(
        "basinmaker/drainage_region_0175_v2-1/finalcat_info_v2-1.zip"
    )
    df = open_shapefile(
        routing_product_shp_path,
    )
    rvh_extractor = BasinMakerExtractor(
        df=df,
        hru_aspect_convention="ArcGIS",
        routing_product_version="2.1",
    )
    rvh_config = rvh_extractor.extract(hru_from_sb=True)
    rvh_config.pop("channel_profile")

    config = BasicRoute(**rvh_config)
    config.write_rv(tmp_path, modelname="routing")
