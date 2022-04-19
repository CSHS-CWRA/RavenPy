from pathlib import Path

import geopandas

from ravenpy.config.commands import SoilProfilesCommand
from ravenpy.extractors.HydrotelShapefileExtractor import HydrotelShapefileExtractor
from ravenpy.extractors.routing_product import RoutingProductShapefileExtractor

pth = "/home/mohammad/Dossier_travail/Hydrotel/DEH/MG24HA/SLSO_MG24HA_2020/physitel/HRU/raven-testdata/famine/hru_Famine_final.zip"
rvh_extractor = RoutingProductShapefileExtractor(
    pth,
    hru_aspect_convention="ArcGIS",
    routing_product_version="2.1",
)

rvh_config = rvh_extractor.extract()


rvh_plus_soil = HydrotelShapefileExtractor(pth, rvh_config)
rvh_plus_soil.extract2(rvh_config)
print(rvh_plus_soil)
