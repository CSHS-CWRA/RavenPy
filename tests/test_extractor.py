from ravenpy.extractors.routing_product import (
    RoutingProductGridWeightExtractor,
    RoutingProductShapefileExtractor,
)
from ravenpy.new_config.emulators.routing import BasicRoute


def test_routing_product_shapefile_extractor(get_file):
    routing_product_shp_path = get_file("raven-routing-sample/finalcat_hru_info.zip")

    rvh_extractor = RoutingProductShapefileExtractor(
        routing_product_shp_path,
        hru_aspect_convention="ArcGIS",
        routing_product_version="1.0",
    )
    rvh_config = rvh_extractor.extract_new_config()
    rvh_config.pop("channel_profiles")

    config = BasicRoute(**rvh_config)
    config.write("/tmp/test_route")
    #      **rvh_config,
    #                 SBGroupPropertyMultiplier=[{"group_name": "Land", "parameter_name":
    # "MANNINGS_N", "mult": 1.0}, {"group_name": "Lakes", "parameter_name": "RESERVOIR_CREST_WIDTH", "mult": 1.0}])
