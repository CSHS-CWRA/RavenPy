from pathlib import Path

import click
import numpy as np

from ravenpy.utilities import gis_import_error_message

try:
    import geopandas as gpd
except (ImportError, ModuleNotFoundError) as e:
    msg = gis_import_error_message.format(Path(__file__).stem)
    raise ImportError(msg) from e


COLUMN_NAMES_CONSTANT_HRU = [
    "SubId",
    "DowSubId",
    "RivSlope",
    "RivLength",
    "BasSlope",
    "BasAspect",
    "BasArea",
    "BkfWidth",
    "BkfDepth",
    "Lake_Cat",
    "HyLakeId",
    "LakeVol",
    "LakeDepth",
    "LakeArea",
    "Laketype",
    "Has_Gauge",
    "MeanElev",
    "FloodP_n",
    "Q_Mean",
    "Ch_n",
    "DrainArea",
    "Strahler",
    "Seg_ID",
    "Seg_order",
    "Max_DEM",
    "Min_DEM",
    "DA_Obs",
    "DA_error",
    "Obs_NM",
    "SRC_obs",
    "centroid_x",
    "centroid_y",
    "HRU_IsLake",
    "Landuse_ID",
    "Soil_ID",
    "Veg_ID",
    "HRU_Area",
    "HRU_ID",
    "LAND_USE_C",
    "VEG_C",
    "SOIL_PROF",
    "HRU_CenX",
    "HRU_CenY",
    "HRU_S_mean",
    "HRU_A_mean",
    "HRU_E_mean",
    "O_ID_1",
    "O_ID_2",
    "geometry",
]

LAND_USE_C_LAKE = "WATER"
LAND_USE_C_LAND = "Landuse_Land_HRU"
SOIL_PROF_LAKE = "LAKE"
SOIL_PROF_LAND = "Soil_Land_HRU"
VEG_C_LAKE = "WATER"
VEG_C_LAND = "Veg_Land_HRU"


@click.command()
@click.argument("input-file", type=click.Path(exists=True))
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    help="Output shapefiles (will create a folder if no extension)",
)
def generate_hrus_from_routing_product(input_file, output):
    """
    Create a new HRU shapefile by splitting every subbasin row of a Routing Product V2.1 shapefile into at least a land HRU and possibly a lake HRU.

    Parameters
    ----------
    input_file : str
        Routing Product V2.1 shapefile (e.g. "drainage_region_0003_v2-1/finalcat_info_v2-1.shp").
    output : str
        Output shapefile (e.g. "hrus.shp").
    """

    def assign_hru_attributes(i_sub, i_hru, hru_id, is_lake_HRU):  # noqa: N803
        # fist copy subbasin attribute to hru table
        for sub_col in subbasin_info.columns:
            if sub_col == "geometry":
                subid = int(subbasin_info["SubId"].values[i_sub])
                # add a  new row in hruinfo
                hru_info.loc[i_hru, sub_col] = np.nan
                # copy geometry to hru
                hru_info.loc[[i_hru], sub_col] = subbasin_info.loc[
                    [subid], sub_col
                ].values
            else:
                hru_info.loc[i_hru, sub_col] = subbasin_info[sub_col].values[i_sub]

        hru_info.loc[i_hru, "HRU_ID"] = hru_id
        hru_info.loc[i_hru, "HRU_CenX"] = subbasin_info["centroid_x"].values[i_sub]
        hru_info.loc[i_hru, "HRU_CenY"] = subbasin_info["centroid_y"].values[i_sub]
        hru_info.loc[i_hru, "HRU_S_mean"] = subbasin_info["BasSlope"].values[i_sub]
        hru_info.loc[i_hru, "HRU_A_mean"] = subbasin_info["BasAspect"].values[i_sub]
        hru_info.loc[i_hru, "HRU_E_mean"] = subbasin_info["MeanElev"].values[i_sub]

        if is_lake_HRU:
            hru_info.loc[i_hru, "HRU_IsLake"] = 1
            hru_info.loc[i_hru, "Soil_ID"] = -1
            hru_info.loc[i_hru, "Veg_ID"] = -1
            hru_info.loc[i_hru, "O_ID_1"] = -1
            hru_info.loc[i_hru, "O_ID_2"] = -1
            hru_info.loc[i_hru, "Landuse_ID"] = -1
            hru_info.loc[i_hru, "HRU_Area"] = subbasin_info["LakeArea"].values[i_sub]
            hru_info.loc[i_hru, "LAND_USE_C"] = LAND_USE_C_LAKE
            hru_info.loc[i_hru, "SOIL_PROF"] = SOIL_PROF_LAKE
            hru_info.loc[i_hru, "VEG_C"] = VEG_C_LAKE
        else:
            hru_info.loc[i_hru, "HRU_IsLake"] = -1
            hru_info.loc[i_hru, "Soil_ID"] = 1
            hru_info.loc[i_hru, "Veg_ID"] = 1
            hru_info.loc[i_hru, "O_ID_1"] = 1
            hru_info.loc[i_hru, "O_ID_2"] = 1
            hru_info.loc[i_hru, "Landuse_ID"] = 1
            hru_info.loc[i_hru, "HRU_Area"] = (
                subbasin_info["BasArea"].values[i_sub]
                - subbasin_info["LakeArea"].values[i_sub]
            )
            hru_info.loc[i_hru, "LAND_USE_C"] = LAND_USE_C_LAND
            hru_info.loc[i_hru, "SOIL_PROF"] = SOIL_PROF_LAND
            hru_info.loc[i_hru, "VEG_C"] = VEG_C_LAND

    if Path(input_file).suffix == ".zip":
        input_file = f"zip://{input_file}"

    subbasin_info = gpd.read_file(input_file)

    if not output:
        p = Path(input_file)
        output_file = p.parent / f"{p.stem}_HRUs.shp"
    else:
        output_file = output

    # create an empty pandas table that each row will be 1 hru, columns will be
    # hru's attributes.
    hru_info = subbasin_info.copy(deep=True)

    # dissolve by subid, the subid will be set as index
    subbasin_info = subbasin_info.dissolve(by="SubId")
    # copy index back to new subid
    subbasin_info["SubId"] = subbasin_info.index

    max_subbsin_id = max(subbasin_info["SubId"].values)

    for col in COLUMN_NAMES_CONSTANT_HRU:
        if col != "geometry":
            hru_info[col] = np.nan

    i_hru = 0

    for i_sub in range(len(subbasin_info)):
        sub_id = subbasin_info["SubId"].values[i_sub]

        if subbasin_info["Lake_Cat"].values[i_sub] > 0:
            # hru id of land will be the corresponding subbasin id
            land_hru_id = sub_id
            assign_hru_attributes(
                i_sub=i_sub,
                i_hru=i_hru,
                hru_id=land_hru_id,
                is_lake_HRU=False,
            )

            i_hru += 1

            # hru id of lake will be subbasin id + max_subbsin_id + 10
            lake_hru_id = land_hru_id + max_subbsin_id + 10
            assign_hru_attributes(
                i_sub=i_sub,
                i_hru=i_hru,
                hru_id=lake_hru_id,
                is_lake_HRU=True,
            )
        else:
            assign_hru_attributes(
                i_sub=i_sub,
                i_hru=i_hru,
                hru_id=sub_id,
                is_lake_HRU=False,
            )

        i_hru += 1

    hru_info = hru_info.loc[hru_info["HRU_ID"] > 0]

    hru_info.to_file(output_file)

    click.echo(f"Created {output_file}")
