from ravenpy.extractors.new_config.routing_product import (
    RoutingProductGridWeightExtractor,
    RoutingProductShapefileExtractor,
)
from ravenpy.new_config.emulators import BasicRoute


def test_routing_product_shapefile_extractor(get_local_testdata, tmp_path):
    routing_product_shp_path = get_local_testdata(
        "raven-routing-sample/finalcat_hru_info.zip"
    )

    rvh_extractor = RoutingProductShapefileExtractor(
        routing_product_shp_path,
        hru_aspect_convention="ArcGIS",
        routing_product_version="1.0",
    )
    rvh_config = rvh_extractor.extract()
    rvh_config.pop("channel_profiles")

    config = BasicRoute(**rvh_config)
    config.write_rv(tmp_path, modelname="routing")

    #      **rvh_config,
    #                 SBGroupPropertyMultiplier=[{"group_name": "Land", "parameter_name":
    # "MANNINGS_N", "mult": 1.0}, {"group_name": "Lakes", "parameter_name": "RESERVOIR_CREST_WIDTH", "mult": 1.0}])
