import warnings
from collections import defaultdict
from pathlib import Path
from typing import List

from ravenpy.utilities import gis_import_error_message

try:
    import geopandas
    from osgeo import __version__ as osgeo_version  # noqa
    from osgeo import ogr, osr  # noqa
except (ImportError, ModuleNotFoundError) as e:
    msg = gis_import_error_message.format(Path(__file__).stem)
    raise ImportError(msg) from e

import netCDF4 as nc4
import numpy as np

from ravenpy.config.commands import (
    ChannelProfileCommand,
    GridWeightsCommand,
    HRUsCommand,
    ReservoirCommand,
    SubBasinsCommand,
)


class RoutingProductShapefileExtractor:

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
    HRU_ASPECT_CONVENTION = "GRASS"  # GRASS | ArcGIS
    ROUTING_PRODUCT_VERSION = "1.0"  # 1.0 | 2.1

    def __init__(
        self,
        shapefile_path,
        hru_aspect_convention=HRU_ASPECT_CONVENTION,
        routing_product_version=ROUTING_PRODUCT_VERSION,
    ):
        if isinstance(shapefile_path, (Path, str)):
            if Path(shapefile_path).suffix == ".zip":
                shapefile_path = f"zip://{shapefile_path}"
            self._df = geopandas.read_file(shapefile_path)
        elif isinstance(shapefile_path, geopandas.GeoDataFrame):
            self._df = shapefile_path

        self.hru_aspect_convention = hru_aspect_convention
        self.routing_product_version = routing_product_version

    def extract(self, model=None):
        """
        This will extract the data from the Routing Product shapefile and
        return it as relevant Raven command data objects.

        Parameters
        ----------
        model : Raven subclass
          Emulator class used to customize HRU objects extracted.

        Returns
        -------
        dict
            "subbasins"
               Sequence of `commands.SubBasinsCommand.Record` objects
            "land_subbasins"
               Sequence of land subbasins ids
            "lake_subbasins"
               Sequence of lake subbasins ids
            "reservoirs"
               Sequence of `commands.ReservoirCommand` objects
            "channel_profiles"
               Sequence of `commands.ChannelProfileCommand` objects
            "hrus"
               Sequence of `commands.HRUsCommand.Record` objects

        """
        self.model_cls = model

        subbasin_recs = []
        land_sb_ids = []
        lake_sb_ids = []
        reservoir_cmds = []  # those are meant to be injected inline in the RVH
        channel_profile_cmds = []  # those are meant to be injected inline in the RVH
        hru_recs = []

        # Collect all subbasin_ids for fast lookup in next loop
        subbasin_ids = {int(row["SubId"]) for _, row in self._df.iterrows()}

        for _, row in self._df.iterrows():

            # HRU
            hru_recs.append(self._extract_hru(row))

            subbasin_id = int(row["SubId"])

            is_lake = False

            lake_field = (
                "IsLake" if self.routing_product_version == "1.0" else "Lake_Cat"
            )

            if row[lake_field] > 0 and row["HRU_IsLake"] > 0:
                lake_sb_ids.append(subbasin_id)
                reservoir_cmds.append(self._extract_reservoir(row))
                is_lake = True
            elif row[lake_field] > 0:
                continue
            else:
                land_sb_ids.append(subbasin_id)

            # Subbasin
            sb = self._extract_subbasin(row, is_lake, subbasin_ids)
            subbasin_recs.append(sb)

            # ChannelProfile
            channel_profile_cmds.append(self._extract_channel_profile(row))

        return dict(
            subbasins=subbasin_recs,
            land_subbasin_ids=land_sb_ids,
            lake_subbasin_ids=lake_sb_ids,
            reservoirs=reservoir_cmds,
            channel_profiles=channel_profile_cmds,
            hrus=hru_recs,
        )

    def _extract_subbasin(self, row, is_lake, subbasin_ids) -> SubBasinsCommand.Record:
        subbasin_id = int(row["SubId"])
        # is_lake = row["HRU_IsLake"] >= 0
        riv_length_field = (
            "Rivlen" if self.routing_product_version == "1.0" else "RivLength"
        )
        river_length_in_kms = 0 if is_lake else round(row[riv_length_field] / 1000, 5)
        # river_slope = max(
        #     row["RivSlope"], RoutingProductShapefileExtractor.MAX_RIVER_SLOPE
        # )
        # downstream_id
        downstream_id = int(row["DowSubId"])
        if downstream_id == subbasin_id:
            downstream_id = -1
        elif downstream_id not in subbasin_ids:
            downstream_id = -1
        has_gauge_field = (
            "IsObs" if self.routing_product_version == "1.0" else "Has_Gauge"
        )
        gauged = row[has_gauge_field] > 0 or (
            is_lake and RoutingProductShapefileExtractor.USE_LAKE_AS_GAUGE
        )
        gauge_id = row["Obs_NM"] if gauged else ""
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

    def _extract_reservoir(self, row) -> ReservoirCommand:
        lake_id = int(row["HyLakeId"])

        return ReservoirCommand(
            subbasin_id=int(row["SubId"]),
            hru_id=int(row["HRU_ID"]),
            name=f"Lake_{lake_id}",
            weir_coefficient=RoutingProductShapefileExtractor.WEIR_COEFFICIENT,
            crest_width=row["BkfWidth"],
            max_depth=row["LakeDepth"],
            lake_area=row["HRU_Area"],
        )

    def _extract_channel_profile(self, row) -> ChannelProfileCommand:
        subbasin_id = int(row["SubId"])
        slope = max(row["RivSlope"], RoutingProductShapefileExtractor.MAX_RIVER_SLOPE)

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

        if RoutingProductShapefileExtractor.USE_MANNING_COEFF:
            mann = channeln
        else:
            mann = RoutingProductShapefileExtractor.MANNING_DEFAULT

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

    def _extract_hru(self, row) -> HRUsCommand.Record:

        aspect = row["HRU_A_mean"]

        if self.hru_aspect_convention == "GRASS":
            aspect -= 360
            if aspect < 0:
                aspect += 360
        elif self.hru_aspect_convention == "ArcGIS":
            aspect = 360 - aspect
        else:
            assert False

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

        if self.model_cls is not None:
            # Instantiate HRUs with emulator specific land_use_class, veg_class and soil_profile names.
            if row["LAND_USE_C"] == "Landuse_Land_HRU":
                return self.model_cls.LandHRU(**attrs)

            if row["LAND_USE_C"] == "Landuse_Lake_HRU":
                return self.model_cls.LakeHRU(**attrs)

        # Instantiate HRUs with generic land_use_class, veg_class and soil_profile names.
        return HRUsCommand.Record(
            land_use_class=row["LAND_USE_C"],
            veg_class=row["VEG_C"],
            soil_profile=row["SOIL_PROF"],
            **attrs,
        )


class RoutingProductGridWeightExtractor:

    """
    The original version of this algorithm can be found at: https://github.com/julemai/GridWeightsGenerator
    """

    DIM_NAMES = ("lon_dim", "lat_dim")
    VAR_NAMES = ("lon", "lat")
    ROUTING_ID_FIELD = "HRU_ID"
    NETCDF_INPUT_FIELD = "NetCDF_col"
    AREA_ERROR_THRESHOLD = 0.05
    CRS_LLDEG = 4326  # EPSG id of lat/lon (deg) coordinate reference system (CRS)
    CRS_CAEA = 3573  # EPSG id of equal-area coordinate reference system (CRS)

    def __init__(
        self,
        input_file_path,
        routing_file_path,
        dim_names=DIM_NAMES,
        var_names=VAR_NAMES,
        routing_id_field=ROUTING_ID_FIELD,
        netcdf_input_field=NETCDF_INPUT_FIELD,
        gauge_ids=None,
        sub_ids=None,
        area_error_threshold=AREA_ERROR_THRESHOLD,
    ):
        self._dim_names = tuple(dim_names)
        self._var_names = tuple(var_names)
        self._routing_id_field = routing_id_field
        self._netcdf_input_field = netcdf_input_field
        self._gauge_ids = gauge_ids or []
        self._sub_ids = sub_ids or []
        self._area_error_threshold = area_error_threshold

        assert not (
            self._gauge_ids and self._sub_ids
        ), "Only one of gauge_ids or sub_ids can be specified"

        self._input_is_netcdf = True

        input_file_path = Path(input_file_path)
        if input_file_path.suffix == ".nc":
            # Note that we cannot use xarray because it complains about variables and dimensions
            # having the same name.
            self._input_data = nc4.Dataset(input_file_path)
        elif input_file_path.suffix == ".zip":
            self._input_data = geopandas.read_file(f"zip://{input_file_path}")
            self._input_is_netcdf = False
        elif input_file_path.suffix == ".shp":
            self._input_data = geopandas.read_file(input_file_path)
            self._input_is_netcdf = False
        else:
            raise ValueError(
                "The input file must be a shapefile (.shp or .zip) or NetCDF"
            )

        routing_file_path = Path(routing_file_path)
        if routing_file_path.suffix == ".zip":
            routing_file_path = f"zip://{routing_file_path}"
        self._routing_data = geopandas.read_file(routing_file_path)

    def extract(self) -> GridWeightsCommand:
        self._prepare_input_data()

        # Read routing data

        # WGS 84 / North Pole LAEA Canada
        self._routing_data = self._routing_data.to_crs(
            epsg=RoutingProductGridWeightExtractor.CRS_CAEA
        )

        def keep_only_valid_downsubid_and_obs_nm(g):
            """
            This function receives a group (g) of routing rows that have been grouped by HRU_ID.
            If there is only one row, there is nothing to do (we return it). If there are more than
            one rows, we want to return the one with the DowSubId != -1 (if there are more than one
            we emit a warning). We also want to search for an Obs_NM value different than -9999 in
            any row of the group, and if we find it we set it in the corresponding field of the
            returning row (if we can't find one we also emit a warning).
            """
            if len(g) == 1:
                return g
            hru_id_field = self._routing_id_field
            row = g[g["DowSubId"] != -1].copy()
            if len(row) > 1:
                row = row[:1].copy()
                warnings.warn(
                    f"More than one row with HRU_ID={row[hru_id_field]} having DowSubId = -1"
                )
            obs_nm = g[g["Obs_NM"] != -9999]
            if not obs_nm.empty:
                row["Obs_NM"] = obs_nm["Obs_NM"].iloc[0]
            else:
                warnings.warn(
                    f"All values of Obs_NM are -9999 for rows with HRU_ID={row[hru_id_field]}"
                )

            return row

        # Remove duplicate HRU_IDs while making sure that we keed relevant DowSubId and Obs_NM values
        self._routing_data = self._routing_data.groupby(self._routing_id_field).apply(
            keep_only_valid_downsubid_and_obs_nm
        )

        # Make sure those are ints
        self._routing_data.SubId = self._routing_data.SubId.astype(int)
        self._routing_data.DowSubId = self._routing_data.DowSubId.astype(int)

        if self._gauge_ids:
            # Extract the SubIDs of the gauges that were specified at input
            self._sub_ids = (
                self._routing_data.loc[self._routing_data.Obs_NM.isin(self._gauge_ids)]
                .SubId.unique()
                .tolist()
            )
            if not self._sub_ids:
                raise ValueError(
                    f"No shapes were found with gauge ID (Obs_NM) in {self._gauge_ids}"
                )

        if self._sub_ids:
            # Here we want to extract the network of connected subbasins by going upstream via their DowSubId,
            # starting from the list supplied by the user (either directly, or via their gauge IDs).. We first
            # build a map of downSubID -> subID for efficient lookup
            downsubid_to_subids = defaultdict(set)
            for _, r in self._routing_data.iterrows():
                downsubid_to_subids[r.DowSubId].add(r.SubId)

            expanded_sub_ids = set(self._sub_ids)
            prev = expanded_sub_ids.copy()
            while True:
                curr = set()
                for did in prev:
                    curr |= downsubid_to_subids[did]
                if curr:
                    expanded_sub_ids |= curr
                    prev = curr
                else:
                    break

            self._sub_ids = expanded_sub_ids

        # Reduce the initial dataset with the target Sub IDs
        if self._sub_ids:
            self._routing_data = self._routing_data[
                self._routing_data.SubId.isin(self._sub_ids)
            ]

        # -------------------------------
        # construct all grid cell polygons
        # -------------------------------

        grid_cell_geom_gpd_wkt = self._compute_grid_cell_polygons()

        # -------------------------------
        # Derive overlay and calculate weights
        # -------------------------------

        grid_weights = []

        for _, row in self._routing_data.iterrows():

            poly = ogr.CreateGeometryFromWkt(row.geometry.to_wkt())

            area_basin = poly.Area()
            # bounding box around basin (for easy check of proximity)
            enve_basin = poly.GetEnvelope()

            area_all = 0.0
            # ncells = 0

            row_grid_weights = []

            for ilat in range(self._nlat):
                for ilon in range(self._nlon):

                    # bounding box around grid-cell (for easy check of proximity)
                    enve_gridcell = grid_cell_geom_gpd_wkt[ilat][ilon].GetEnvelope()

                    grid_is_close = self._check_proximity_of_envelops(
                        enve_gridcell, enve_basin
                    )

                    # this check decreases runtime DRASTICALLY (from ~6h to ~1min)
                    if not grid_is_close:
                        continue

                    # grid_cell_area = grid_cell_geom_gpd_wkt[ilat][ilon].Area()

                    # "fake" buffer to avoid invalid polygons and weirdos dumped by ArcGIS
                    inter = grid_cell_geom_gpd_wkt[ilat][ilon].Intersection(
                        poly.Buffer(0.0)
                    )

                    area_intersect = inter.Area()

                    area_all += area_intersect

                    if area_intersect > 0:
                        hru_id = int(row[self._routing_id_field])
                        cell_id = ilat * self._nlon + ilon
                        weight = area_intersect / area_basin
                        row_grid_weights.append((hru_id, cell_id, weight))

            # mismatch between area of subbasin (routing product) and sum of all contributions of grid cells (model output)
            error = (area_basin - area_all) / area_basin

            if abs(error) > self._area_error_threshold and area_basin > 500000.0:
                # record all basins with errors larger 5% (if basin is larger than 0.5 km2)
                # error_dict[int(ibasin[key_colname])] = [error, area_basin]
                grid_weights += row_grid_weights

            else:
                # adjust such that weights sum up to 1.0
                for hru_id, cell_id, weight in row_grid_weights:
                    corrected_weight = weight * 1.0 / (1.0 - error)
                    grid_weights.append((hru_id, cell_id, corrected_weight))

                # if error < 1.0:
                #     area_all *= 1.0 / (1.0 - error)
                # error = 0.0

        return GridWeightsCommand(
            number_hrus=len(self._routing_data),
            number_grid_cells=self._nlon * self._nlat,
            data=tuple(grid_weights),
        )

    def _prepare_input_data(self):

        if self._input_is_netcdf:

            # Raven numbering is:
            #
            #      [      1      2      3   ...     1*nlon
            #        nlon+1 nlon+2 nlon+3   ...     2*nlon
            #           ...    ...    ...   ...     ...
            #           ...    ...    ...   ...  nlat*nlon ]
            #
            # --> Making sure shape of lat/lon fields is like that
            #

            lon_var = self._input_data.variables[self._var_names[0]]
            lon_dims = lon_var.dimensions
            lon_var = lon_var[:]

            lat_var = self._input_data.variables[self._var_names[1]]
            lat_dims = lat_var.dimensions
            lat_var = lat_var[:]

            if len(lon_dims) == 2 and len(lat_dims) == 2:
                if lon_dims == self._dim_names:
                    lon_var = np.transpose(lon_var)

                if lat_dims == self._dim_names:
                    lat_var = np.transpose(lat_var)

            elif len(lon_dims) == 1 and len(lat_dims) == 1:
                # Coord vars are 1d, make them 2d
                lon_var_size = lon_var.size
                lon_var = np.tile(lon_var, (lat_var.size, 1))
                lat_var = np.transpose(np.tile(lat_var, (lon_var_size, 1)))
                assert lon_var.shape == lat_var.shape

            else:
                raise ValueError(
                    "The coord variables of the input data must have the same number of dimensions (either 1 or 2)"
                )

            self._lath, self._lonh = self._create_gridcells_from_centers(
                lat_var, lon_var
            )

            self._nlon = np.shape(lon_var)[1]
            self._nlat = np.shape(lat_var)[0]

        else:

            # input data is a shapefile

            self._input_data = self._input_data.to_crs(
                epsg=RoutingProductGridWeightExtractor.CRS_CAEA
            )

            self._nlon = 1  # only for consistency

            # number of shapes in model "discretization" shapefile (not routing toolbox shapefile)
            self._nlat = self._input_data.geometry.count()  # only for consistency

    def _compute_grid_cell_polygons(self):

        grid_cell_geom_gpd_wkt: List[List[List[ogr.Geometry]]] = [
            [[] for ilon in range(self._nlon)] for ilat in range(self._nlat)
        ]

        if self._input_is_netcdf:

            lath = self._lath
            lonh = self._lonh

            for ilat in range(self._nlat):

                for ilon in range(self._nlon):

                    # -------------------------
                    # EPSG:3035   needs a swap before and after transform ...
                    # -------------------------
                    # gridcell_edges = [ [lath[ilat,ilon]    , lonh[ilat,  ilon]    ],            # for some reason need to switch lat/lon that transform works
                    #                    [lath[ilat+1,ilon]  , lonh[ilat+1,ilon]    ],
                    #                    [lath[ilat+1,ilon+1], lonh[ilat+1,ilon+1]  ],
                    #                    [lath[ilat,ilon+1]  , lonh[ilat,  ilon+1]  ]]

                    # tmp = self._routing_data_to_geometry(gridcell_edges, epsg=crs_caea)
                    # tmp.SwapXY()              # switch lat/lon back
                    # grid_cell_geom_gpd_wkt[ilat][ilon] = tmp

                    # -------------------------
                    # EPSG:3573   does not need a swap after transform ... and is much faster than transform with EPSG:3035
                    # -------------------------
                    #
                    # Windows            Python 3.8.5 GDAL 3.1.3 --> lat/lon (Ming)
                    # MacOS 10.15.6      Python 3.8.5 GDAL 3.1.3 --> lat/lon (Julie)
                    # Graham             Python 3.8.2 GDAL 3.0.4 --> lat/lon (Julie)
                    # Graham             Python 3.6.3 GDAL 2.2.1 --> lon/lat (Julie)
                    # Ubuntu 18.04.2 LTS Python 3.6.8 GDAL 2.2.3 --> lon/lat (Etienne)
                    #
                    if osgeo_version < "3.0":
                        gridcell_edges = [
                            [
                                lonh[ilat, ilon],
                                lath[ilat, ilon],
                            ],  # for some reason need to switch lat/lon that transform works
                            [lonh[ilat + 1, ilon], lath[ilat + 1, ilon]],
                            [lonh[ilat + 1, ilon + 1], lath[ilat + 1, ilon + 1]],
                            [lonh[ilat, ilon + 1], lath[ilat, ilon + 1]],
                        ]
                    else:
                        gridcell_edges = [
                            [
                                lath[ilat, ilon],
                                lonh[ilat, ilon],
                            ],  # for some reason lat/lon order works
                            [lath[ilat + 1, ilon], lonh[ilat + 1, ilon]],
                            [lath[ilat + 1, ilon + 1], lonh[ilat + 1, ilon + 1]],
                            [lath[ilat, ilon + 1], lonh[ilat, ilon + 1]],
                        ]

                    tmp = self._shape_to_geometry(
                        gridcell_edges, epsg=RoutingProductGridWeightExtractor.CRS_CAEA
                    )
                    grid_cell_geom_gpd_wkt[ilat][ilon] = tmp

        else:

            for ishape in range(self._nlat):

                idx = np.where(self._input_data[self._netcdf_input_field] == ishape)[0]
                if len(idx) == 0:
                    # print(
                    #     "Polygon ID = {} not found in '{}'. Numbering of shapefile attribute '{}' needs to be [0 ... {}-1].".format(
                    #         ishape, input_file, key_colname_model, nshapes
                    #     )
                    # )
                    raise ValueError("Polygon ID not found.")
                if len(idx) > 1:
                    # print(
                    #     "Polygon ID = {} found multiple times in '{}' but needs to be unique. Numbering of shapefile attribute '{}' needs to be [0 ... {}-1].".format(
                    #         ishape, input_file, key_colname_model, nshapes
                    #     )
                    # )
                    raise ValueError("Polygon ID not unique.")
                idx = idx[0]
                poly = self._input_data.loc[idx].geometry
                grid_cell_geom_gpd_wkt[ishape][0] = ogr.CreateGeometryFromWkt(
                    poly.to_wkt()
                ).Buffer(
                    0.0
                )  # We add an empty buffer here to fix problems with bad polygon topology (actually caused by ESRI's historical incompetence)

        return grid_cell_geom_gpd_wkt

    def _create_gridcells_from_centers(self, lat, lon):

        # create array of edges where (x,y) are always center cells
        nlon = np.shape(lon)[1]
        nlat = np.shape(lat)[0]
        lonh = np.empty((nlat + 1, nlon + 1), dtype=np.float64)
        lath = np.empty((nlat + 1, nlon + 1), dtype=np.float64)
        tmp1 = [
            [(lat[ii + 1, jj + 1] - lat[ii, jj]) / 2 for jj in range(nlon - 1)]
            + [(lat[ii + 1, nlon - 1] - lat[ii, nlon - 2]) / 2]
            for ii in range(nlat - 1)
        ]
        tmp2 = [
            [(lon[ii + 1, jj + 1] - lon[ii, jj]) / 2 for jj in range(nlon - 1)]
            + [(lon[ii + 1, nlon - 1] - lon[ii, nlon - 2]) / 2]
            for ii in range(nlat - 1)
        ]
        dlat = np.array(tmp1 + [tmp1[-1]])
        dlon = np.array(tmp2 + [tmp2[-1]])
        lonh[0:nlat, 0:nlon] = lon - dlon
        lath[0:nlat, 0:nlon] = lat - dlat

        # make lat and lon one column and row wider such that all
        lonh[nlat, 0:nlon] = lonh[nlat - 1, 0:nlon] + (
            lonh[nlat - 1, 0:nlon] - lonh[nlat - 2, 0:nlon]
        )
        lath[nlat, 0:nlon] = lath[nlat - 1, 0:nlon] + (
            lath[nlat - 1, 0:nlon] - lath[nlat - 2, 0:nlon]
        )
        lonh[0:nlat, nlon] = lonh[0:nlat, nlon - 1] + (
            lonh[0:nlat, nlon - 1] - lonh[0:nlat, nlon - 2]
        )
        lath[0:nlat, nlon] = lath[0:nlat, nlon - 1] + (
            lath[0:nlat, nlon - 1] - lath[0:nlat, nlon - 2]
        )
        lonh[nlat, nlon] = lonh[nlat - 1, nlon - 1] + (
            lonh[nlat - 1, nlon - 1] - lonh[nlat - 2, nlon - 2]
        )
        lath[nlat, nlon] = lath[nlat - 1, nlon - 1] + (
            lath[nlat - 1, nlon - 1] - lath[nlat - 2, nlon - 2]
        )

        return [lath, lonh]

    def _shape_to_geometry(self, shape_from_jsonfile, epsg=None):

        # converts shape read from shapefile to geometry
        # epsg :: integer EPSG code

        ring_shape = ogr.Geometry(ogr.wkbLinearRing)

        for ii in shape_from_jsonfile:
            ring_shape.AddPoint_2D(ii[0], ii[1])
        # close ring
        ring_shape.AddPoint_2D(shape_from_jsonfile[0][0], shape_from_jsonfile[0][1])

        poly_shape = ogr.Geometry(ogr.wkbPolygon)
        poly_shape.AddGeometry(ring_shape)

        if epsg:
            source = osr.SpatialReference()
            # usual lat/lon projection
            source.ImportFromEPSG(RoutingProductGridWeightExtractor.CRS_LLDEG)

            target = osr.SpatialReference()
            target.ImportFromEPSG(epsg)  # any projection to convert to

            transform = osr.CoordinateTransformation(source, target)
            poly_shape.Transform(transform)

        return poly_shape

    def _check_proximity_of_envelops(self, gridcell_envelop, shape_envelop):

        # checks if two envelops are in proximity (intersect)

        # minX  --> env[0]
        # maxX  --> env[1]
        # minY  --> env[2]
        # maxY  --> env[3]

        return (
            gridcell_envelop[0] <= shape_envelop[1]
            and gridcell_envelop[1] >= shape_envelop[0]
            and gridcell_envelop[2] <= shape_envelop[3]
            and gridcell_envelop[3] >= shape_envelop[2]
        )

    def _check_gridcell_in_proximity_of_shape(
        self, gridcell_edges, shape_from_jsonfile
    ):

        # checks if a grid cell falls into the bounding box of the shape
        # does not mean it intersects but it is a quick and cheap way to
        # determine cells that might intersect

        # gridcell_edges = [(lon1,lat1),(lon2,lat2),(lon3,lat3),(lon4,lat4)]
        # shape_from_jsonfile

        min_lat_cell = np.min([ii[1] for ii in gridcell_edges])
        max_lat_cell = np.max([ii[1] for ii in gridcell_edges])
        min_lon_cell = np.min([ii[0] for ii in gridcell_edges])
        max_lon_cell = np.max([ii[0] for ii in gridcell_edges])

        lat_shape = np.array(
            [icoord[1] for icoord in shape_from_jsonfile]
        )  # is it lat???
        lon_shape = np.array(
            [icoord[0] for icoord in shape_from_jsonfile]
        )  # is it lon???

        min_lat_shape = np.min(lat_shape)
        max_lat_shape = np.max(lat_shape)
        min_lon_shape = np.min(lon_shape)
        max_lon_shape = np.max(lon_shape)

        return (
            min_lat_cell <= max_lat_shape
            and max_lat_cell >= min_lat_shape
            and min_lon_cell <= max_lon_shape
            and max_lon_cell >= min_lon_shape
        )
