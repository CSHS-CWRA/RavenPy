import warnings
from collections import defaultdict
from pathlib import Path
from typing import List

import geopandas

from ravenpy.config.commands import SoilProfilesCommand


class HydrotelShapefileExtractor:
    def __init__(self, shapefile_path, routingproductoutput):
        self.routingproductoutput = routingproductoutput
        if isinstance(shapefile_path, (Path, str)):
            if Path(shapefile_path).suffix == ".zip":
                shapefile_path = f"zip://{shapefile_path}"
            self._df = geopandas.read_file(shapefile_path)
        elif isinstance(shapefile_path, geopandas.GeoDataFrame):
            self._df = shapefile_path

    def extract2(self, routingproductoutput) -> SoilProfilesCommand.Record:
        # This method extracts the soil profile data (layer name, thickness) from the shapefile.

        soil_recs = set()
        for _, row in self._df.iterrows():
            soil_recs.add(
                self._extract_soillayer_thickness(row)
            )  # here we find the unique rows
        soils = list(soil_recs)
        soilinfo = []
        keyss = ["soilprofile"]
        for soil in soils:
            soilinfo.append(
                SoilProfilesCommand.Record(
                    soil[0], ["TOPSOIL", "FAST_RES", "SLOW_RES"], soil[2:]
                )
            )
        ss = dict.fromkeys(keyss, soilinfo)
        return routingproductoutput.update(ss)

    def _extract_soillayer_thickness(self, row):

        items = [
            "TOPSOIL",
            "Nhorizons",
            "th1",
            "th2",
            "th3",
        ]  # all 3 layers have the same soil type (Hydrotel configurations)
        return tuple(row[item] for item in items)
