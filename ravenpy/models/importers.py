import geopandas

from ravenpy.models.state import HRURecord, SubbasinLakeRecord, SubbasinRecord


class RoutingProductShapefileImporter:

    """
    This is a class to encapsulate the logic of converting the Routing
    Product into the required data structures to generate the RVH file
    format.

    """

    MAX_RIVER_SLOPE = 0.00001
    USE_LAKE_AS_GAUGE = False
    WEIR_COEFFICIENT = 0.6
    USE_LAND_AS_GAUGE = False

    def __init__(self, path):
        self._df = geopandas.read_file(path)

    def extract(self):
        sbs, groups, lakes = self._extract_subbasins_and_lakes()
        hrus = self._extract_hrus()

        return sbs, groups, lakes, hrus

    def _extract_lake(self, row):
        lake_id = int(row["HyLakeId"])

        return SubbasinLakeRecord(
            subbasin_id=int(row["SubId"]),
            hru_id=int(row["HRU_ID"]),
            name=f"Lake_{lake_id}",
            weir_coefficient=RoutingProductShapefileImporter.WEIR_COEFFICIENT,
            crest_width=row["BkfWidth"],
            max_depth=row["LakeDepth"],
            lake_area=row["HRU_Area"],
        )

    def _extract_subbasins_and_lakes(self):
        # Collect all subbasin_ids for fast lookup in next loop
        subbasin_ids = {int(row["SubId"]) for _, row in self._df.iterrows()}

        subbasins = []
        subbasin_groups = {
            "land": [],
            "lake": [],
        }
        lakes = []

        # Here we only consider records with unique SubIds
        for _, row in self._df.drop_duplicates("SubId", keep="first").iterrows():
            subbasin_id = int(row["SubId"])

            if row["IsLake"] >= 0:
                river_length_in_kms = 0
                subbasin_groups["lake"].append(subbasin_id)
                lakes.append(self._extract_lake(row))
            else:
                river_length_in_kms = row["Rivlen"] / 1000
                subbasin_groups["land"].append(subbasin_id)

            river_slope = max(
                row["RivSlope"], RoutingProductShapefileImporter.MAX_RIVER_SLOPE
            )

            # downstream_id
            downstream_id = int(row["DowSubId"])
            if downstream_id == subbasin_id:
                downstream_id = -1
            elif downstream_id not in subbasin_ids:
                downstream_id = -1

            gauged = row["IsObs"] > 0 or (
                row["IsLake"] >= 0 and RoutingProductShapefileImporter.USE_LAKE_AS_GAUGE
            )

            rec = SubbasinRecord(
                subbasin_id=subbasin_id,
                name=f"sub{subbasin_id}",
                downstream_id=downstream_id,
                profile=f"chn_{subbasin_id}",
                reach_length=river_length_in_kms,
                gauged=gauged,
            )
            subbasins.append(rec)

        return subbasins, subbasin_groups, lakes

    def _extract_hrus(self):
        hrus = []
        for _, row in self._df.iterrows():
            hru = HRURecord(
                hru_id=int(row["HRU_ID"]),
                area=row["HRU_Area"],
                elevation=row["HRU_E_mean"],
                latitude=row["HRU_CenY"],
                longitude=row["HRU_CenX"],
                subbasin_id=int(row["SubId"]),
                land_use_class=row["LAND_USE_C"],
                veg_class=row["VEG_C"],
                soil_profile=row["SOIL_PROF"],
                aquifer_profile="[NONE]",
                terrain_class="[NONE]",
                slope=row["HRU_S_mean"],
                aspect=row["HRU_A_mean"],
            )
            hrus.append(hru)

        return hrus
