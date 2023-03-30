from ravenpy.extractors.new_config.routing_product import (
    BasinMakerExtractor,
    GridWeightExtractor,
    open_shapefile,
)
from ravenpy.new_config.emulators import BasicRoute


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
    rvh_config = rvh_extractor.extract()
    rvh_config.pop("channel_profiles")

    config = BasicRoute(**rvh_config)
    config.write_rv(tmp_path, modelname="routing")

    #      **rvh_config,
    #                 SBGroupPropertyMultiplier=[{"group_name": "Land", "parameter_name":
    # "MANNINGS_N", "mult": 1.0}, {"group_name": "Lakes", "parameter_name": "RESERVOIR_CREST_WIDTH", "mult": 1.0}])
