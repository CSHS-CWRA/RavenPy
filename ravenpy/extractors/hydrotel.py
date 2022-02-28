from pathlib import Path

import geopandas

from ravenpy.config.commands import (
    ChannelProfileCommand,
    HRUsCommand,
    ReservoirCommand,
    SubBasinsCommand,
)

MAX_RIVER_SLOPE = 0.00001
USE_LAKE_AS_GAUGE = False
WEIR_COEFFICIENT = 0.6
USE_MANNING_COEFF = False
MANNING_DEFAULT = 0.035
HRU_ASPECT_CONVENTION = "GRASS"  # GRASS | ArcGIS


def extract_hru(row, aspect_convention="GRASS") -> HRUsCommand.Record:

    aspect = row["HRU_A_mean"]

    if HRU_ASPECT_CONVENTION == "GRASS":
        aspect -= 360
        if aspect < 0:
            aspect += 360
    else:
        aspect = 360 - aspect

    attrs = dict(
        hru_id=int(row["HRU_ID"]),
        area=row["HRU_Area"] / 1_000_000,
        elevation=row["HRU_E_mean"],
        latitude=row["HRU_CenY"],
        longitude=row["HRU_CenX"],
        subbasin_id=int(row["SubId"]),
        aquifer_profile="[NONE]",
        terrain_class="[NONE]",
        slope=row["HRU_S_mean"],
        aspect=aspect,
    )

    # Instantiate HRUs with generic land_use_class, veg_class and soil_profile names.
    return HRUsCommand.Record(
        land_use_class=row["LAND_USE_C"],
        veg_class=row["VEG_C"],
        soil_profile=row["SOIL_PROF"],
        **attrs,
    )


def extract_subbasin(row, subbasin_ids) -> SubBasinsCommand.Record:
    subbasin_id = int(row["SubId"])
    # river_length_in_kms = 0 if is_lake else round(row["RivLength"] / 1000, 5)
    river_length_in_kms = round(row["RivLength"] / 1000, 5)
    downstream_id = int(row["DowSubId"])
    if downstream_id == subbasin_id:
        downstream_id = -1
    elif downstream_id not in subbasin_ids:
        downstream_id = -1
    # has_gauge_field = "IsObs" if self.routing_product_version == "1.0" else "Has_Gauge"
    # gauged = row[has_gauge_field] > 0 or (
    #     is_lake and RoutingProductShapefileExtractor.USE_LAKE_AS_GAUGE
    # )
    # gauge_id = row["Obs_NM"] if gauged else ""
    gauged = False
    gauge_id = ""
    rec = SubBasinsCommand.Record(
        subbasin_id=subbasin_id,
        name=f"sub_{subbasin_id}",
        downstream_id=downstream_id,
        profile=f"chn_{subbasin_id}",
        reach_length=river_length_in_kms,
        gauged=gauged,
        gauge_id=gauge_id,
    )

    return rec


def extract_reservoir(row) -> ReservoirCommand:
    lake_id = int(row["HyLakeId"])

    return ReservoirCommand(
        subbasin_id=int(row["SubId"]),
        hru_id=int(row["HRU_ID"]),
        name=f"Lake_{lake_id}",
        weir_coefficient=WEIR_COEFFICIENT,
        crest_width=row["BkfWidth"],
        max_depth=row["LakeDepth"],
        lake_area=row["HRU_Area"],
    )


def extract_channel_profile(row) -> ChannelProfileCommand:
    subbasin_id = int(row["SubId"])
    slope = max(row["RivSlope"], MAX_RIVER_SLOPE)

    # SWAT: top width of channel when filled with water; bankfull width W_bnkfull
    channel_width = max(row["BkfWidth"], 1)
    # SWAT: depth of water in channel when filled to top of bank
    channel_depth = max(row["BkfDepth"], 1)
    channel_elev = row["MeanElev"]
    floodn = row["FloodP_n"]
    channeln = row["Ch_n"]

    # channel profile calculations are based on theory SWAT model is based on
    # see: https://swat.tamu.edu/media/99192/swat2009-theory.pdf
    #      --> "Channel Characteristics" p. 429 ff

    # inverse of channel side slope; channel sides assumed to have 2:1 run to rise ratio
    zch = 2
    # river side width
    sidwd = zch * channel_depth
    # river bottom width W_btm
    botwd = channel_width - 2 * sidwd

    # if derived bottom width is negative, set bottom width to 0.5*bankfull width and recalculate zch
    if botwd < 0:
        botwd = 0.5 * channel_width
        sidwd = 0.5 * 0.5 * channel_width
        zch = (channel_width - botwd) / (2 * channel_depth)

    # inverse of floodplain side slope; flood plain side slopes assumed to have 4:1 run to rise ratio
    zfld = 4 + channel_elev
    # floodplain bottom width
    zbot = channel_elev - channel_depth
    # floodplain side width
    sidwdfp = 4 / 0.25

    # geometry of the channel and floodplain
    # (see figure 7:1-2 in SWAT theory document)
    survey_points = (
        (0, zfld),
        (sidwdfp, channel_elev),
        (sidwdfp + 2 * channel_width, channel_elev),
        (sidwdfp + 2 * channel_width + sidwd, zbot),
        (sidwdfp + 2 * channel_width + sidwd + botwd, zbot),
        (sidwdfp + 2 * channel_width + 2 * sidwd + botwd, channel_elev),
        (sidwdfp + 4 * channel_width + 2 * sidwd + botwd, channel_elev),
        (2 * sidwdfp + 4 * channel_width + 2 * sidwd + botwd, zfld),
    )

    if USE_MANNING_COEFF:
        mann = channeln
    else:
        mann = MANNING_DEFAULT

    # roughness zones of channel and floodplain
    roughness_zones = (
        (0, floodn),
        (sidwdfp + 2 * channel_width, mann),
        (sidwdfp + 2 * channel_width + 2 * sidwd + botwd, floodn),
    )

    return ChannelProfileCommand(
        name=f"chn_{subbasin_id}",
        bed_slope=slope,
        survey_points=survey_points,
        roughness_zones=roughness_zones,
    )


def extract_rv_objects_from_shapefile(shapefile_path):
    if Path(shapefile_path).suffix == ".zip":
        shapefile_path = f"zip://{shapefile_path}"
    df = geopandas.read_file(shapefile_path)

    subbasin_recs = {}
    hru_recs = []
    reservoir_cmds = []
    channel_profile_cmds = []

    # Collect all subbasin_ids for fast lookup in extract_subbasin
    subbasin_ids = {int(row["SubId"]) for _, row in df.iterrows()}

    for _, row in df.iterrows():

        # Each row corresponds to an HRU
        hru_recs.append(extract_hru(row))

        sb = extract_subbasin(row, subbasin_ids)
        subbasin_recs[sb.subbasin_id] = sb

        if row["Lake_Cat"] > 0 and row["HRU_IsLake"] > 0:
            reservoir_cmds.append(extract_reservoir(row))

        channel_profile_cmds.append(extract_channel_profile(row))

    return dict(
        subbasins=list(subbasin_recs.values()),
        reservoirs=reservoir_cmds,
        channel_profiles=channel_profile_cmds,
        hrus=hru_recs,
    )
