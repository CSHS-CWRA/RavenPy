import geopandas

from .commands import (
    ChannelProfileCommand,
    HRUsCommand,
    HRUsCommandRecord,
    ReservoirCommand,
    SubBasinGroupCommand,
    SubBasinsCommand,
    SubBasinsCommandRecord,
)


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
    USE_MANNING_COEFF = False
    MANNING_DEFAULT = 0.035

    def __init__(self, path):
        self._df = geopandas.read_file(path)

    def extract(self):
        """
        This will extract the data from the Routing Product shapefile and
        return it as relevant commands.

        Returns
        -------
        `commands.SubBasinsCommand`
        `commands.SubBasinGroup`
        `commands.SubBasinGroup`
        list of `commands.ReservoirCommand`
        list of `commands.ChannelProfileCommand`
        `commands.HRUsCommand`

        """

        subbasin_recs = []
        land_sb_ids = []
        lake_sb_ids = []
        reservoir_cmds = []  # those are meant to be injected inline in the RVH
        channel_profile_cmds = []  # those are meant to be injected inline in the RVH
        hru_recs = []

        # Collect all subbasin_ids for fast lookup in next loop
        subbasin_ids = {int(row["SubId"]) for _, row in self._df.iterrows()}
        subbasin_id_accum = set()

        for _, row in self._df.iterrows():

            # HRU
            hru_recs.append(self._extract_hru(row))

            subbasin_id = int(row["SubId"])

            # We only want to process the first row with a given SubId (we ignore the other ones)
            if subbasin_id in subbasin_id_accum:
                continue
            subbasin_id_accum.add(subbasin_id)

            # Subbasin
            sb, is_lake = self._extract_subbasin(row, subbasin_ids)
            subbasin_recs.append(sb)

            if is_lake:
                lake_sb_ids.append(subbasin_id)
                reservoir_cmds.append(self._extract_reservoir(row))
            else:
                land_sb_ids.append(subbasin_id)

            # ChannelProfile
            channel_profile_cmds.append(self._extract_channel_profile(row))

        return (
            SubBasinsCommand(subbasin_recs),
            SubBasinGroupCommand("land", land_sb_ids),
            SubBasinGroupCommand("lake", lake_sb_ids),
            reservoir_cmds,
            channel_profile_cmds,
            HRUsCommand(hru_recs),
        )

    def _extract_subbasin(self, row, subbasin_ids) -> SubBasinsCommandRecord:
        subbasin_id = int(row["SubId"])
        is_lake = row["IsLake"] >= 0
        river_length_in_kms = 0 if is_lake else row["Rivlen"] / 1000
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
        rec = SubBasinsCommandRecord(
            subbasin_id=subbasin_id,
            name=f"sub_{subbasin_id}",
            downstream_id=downstream_id,
            profile=f"chn_{subbasin_id}",
            reach_length=river_length_in_kms,
            gauged=gauged,
        )

        return rec, is_lake

    def _extract_reservoir(self, row) -> ReservoirCommand:
        lake_id = int(row["HyLakeId"])

        return ReservoirCommand(
            subbasin_id=int(row["SubId"]),
            hru_id=int(row["HRU_ID"]),
            name=f"Lake_{lake_id}",
            weir_coefficient=RoutingProductShapefileImporter.WEIR_COEFFICIENT,
            crest_width=row["BkfWidth"],
            max_depth=row["LakeDepth"],
            lake_area=row["HRU_Area"],
        )

    def _extract_channel_profile(self, row) -> ChannelProfileCommand:
        subbasin_id = int(row["SubId"])
        slope = max(row["RivSlope"], RoutingProductShapefileImporter.MAX_RIVER_SLOPE)

        channel_width = row["BkfWidth"]
        channel_depth = row["BkfDepth"]
        channel_elev = row["MeanElev"]
        floodn = row["FloodP_n"]
        channeln = row["Ch_n"]

        zch = 2
        sidwd = zch * channel_depth  # river side width
        botwd = channel_width - 2 * sidwd  # river
        if botwd < 0:
            botwd = 0.5 * channel_width
            sidwd = 0.5 * 0.5 * channel_width
            zch = (channel_width - botwd) / 2 / channel_depth

        zfld = 4 + channel_elev
        zbot = channel_elev - channel_depth
        sidwdfp = 4 / 0.25

        survey_points = [
            (0, zfld),
            (sidwdfp, channel_elev),
            (sidwdfp + 2 * channel_width, channel_elev),
            (sidwdfp + 2 * channel_width + sidwd, zbot),
            (sidwdfp + 2 * channel_width + sidwd + botwd, zbot),
            (sidwdfp + 2 * channel_width + 2 * sidwd + botwd, channel_elev),
            (sidwdfp + 4 * channel_width + 2 * sidwd + botwd, channel_elev),
            (2 * sidwdfp + 4 * channel_width + 2 * sidwd + botwd, zfld),
        ]

        if RoutingProductShapefileImporter.USE_MANNING_COEFF:
            mann = channeln
        else:
            mann = RoutingProductShapefileImporter.MANNING_DEFAULT

        roughness_zones = [
            (0, floodn),
            (sidwdfp + 2 * channel_width, mann),
            (sidwdfp + 2 * channel_width + 2 * sidwd + botwd, floodn),
        ]

        return ChannelProfileCommand(
            name=f"chn_{subbasin_id}",
            bed_slope=slope,
            survey_points=survey_points,
            roughness_zones=roughness_zones,
        )

    def _extract_hru(self, row) -> HRUsCommandRecord:
        return HRUsCommandRecord(
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
