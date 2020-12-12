import geopandas
import netCDF4 as nc4
import numpy as np
import xarray

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


class GridWeights:
    def __init__(self):
        # nc_in = xarray.open_dataset(
        #     "/home/christian/ouranos/raven/doc/Lievre/input/VIC_streaminputs.nc"
        # )
        dimname = ["lon", "lat"]
        varname = ["lon", "lat"]

        routinginfo = "/home/christian/ouranos/raven/doc/Lievre/maps/LievreHRUs.shp"

        crs_lldeg = 4326  # EPSG id of lat/lon (deg) coordinate referenence system (CRS)
        crs_caea = 3573  # EPSG id of equal-area    coordinate referenence system (CRS)

        key_colname = "HRU_ID"

        doall = False

        basin = "1802"

        nc_in = nc4.Dataset(
            "/home/christian/ouranos/raven/doc/Lievre/input/VIC_streaminputs.nc"
        )

        # print(nc_in)

        lon = nc_in.variables[varname[0]][:]
        lon_dims = nc_in.variables[varname[0]].dimensions
        lat = nc_in.variables[varname[1]][:]
        lat_dims = nc_in.variables[varname[1]].dimensions
        nc_in.close()

        # Raven numbering is:
        #
        #      [      1      2      3   ...     1*nlon
        #        nlon+1 nlon+2 nlon+3   ...     2*nlon
        #           ...    ...    ...   ...     ...
        #           ...    ...    ...   ...  nlat*nlon ]
        #
        # --> Making sure shape of lat/lon fields is like that
        #
        if np.all(np.array(lon_dims) == dimname[::1]):
            lon = np.transpose(lon)
            print(
                '   >>> switched order of dimensions for variable "{0}"'.format(
                    varname[0]
                )
            )
        elif np.all(np.array(lon_dims) == dimname[::-1]):
            print(
                '   >>> order of dimensions correct for variable "{0}"'.format(
                    varname[0]
                )
            )
        else:
            print(
                "   >>> Dimensions found {0} does not match the dimension names specified with (-d): {1}".format(
                    lon_dims, dimname
                )
            )
            raise ValueError("STOP")

        if np.all(np.array(lat_dims) == dimname[::1]):
            lat = np.transpose(lat)
            print(
                '   >>> switched order of dimensions for variable "{0}"'.format(
                    varname[1]
                )
            )
        elif np.all(np.array(lat_dims) == dimname[::-1]):
            print(
                '   >>> order of dimensions correct for variable "{0}"'.format(
                    varname[1]
                )
            )
        else:
            print(
                "   >>> Dimensions found {0} does not match the dimension names specified with (-d): {1}".format(
                    lat_dims, dimname
                )
            )
            raise ValueError("STOP")

        lath, lonh = self.create_gridcells_from_centers(lat, lon)

        nlon = np.shape(lon)[1]
        nlat = np.shape(lat)[0]

        # -------------------------------
        # Read Basin shapes and all subbasin-shapes
        # -------------------------------
        print(" ")
        print("   (2) Reading routing toolbox data ...")

        shape = geopandas.read_file(routinginfo)
        # shape     = shape.to_crs(epsg=crs_lldeg)        # this is lat/lon in degree
        shape = shape.to_crs(epsg=crs_caea)  # WGS 84 / North Pole LAEA Canada

        # check that key column contains only unique values
        keys = np.array(list(shape[key_colname]))
        # keys_uniq = np.unique(keys)
        # if len(keys_uniq) != len(keys):
        #     raise ValueError("The attribute of the shapefile set to contain only unique identifiers ('{}') does contain duplicate keys. Please specify another column (option -c '<col_name>') and use the option to process all records contained in the shapefile (-a).".format(key_colname))

        # select only relevant basins/sub-basins
        if not (doall):

            if not (basin is None):  # if gauge ID is given

                basins = [bb.strip() for bb in basin.split(",")]
                idx_basins = [list(np.where(shape["Obs_NM"] == bb)[0]) for bb in basins]

                # find corresponding SubId
                SubId = [np.int(shape.loc[idx_basin].SubId) for idx_basin in idx_basins]
                print("   >>> found gauge at SubId = ", SubId)

            if not (SubId is None):  # if routing toolbox basin ID is given

                old_SubIds = []
                for SI in SubId:

                    old_SubId = []
                    new_SubId = [SI]

                    while len(new_SubId) > 0:

                        old_SubId.append(new_SubId)
                        new_SubId = [
                            list(
                                shape.loc[(np.where(shape["DowSubId"] == ii))[0]].SubId
                            )
                            for ii in new_SubId
                        ]  # find all upstream catchments of these new basins
                        new_SubId = list(
                            np.unique(
                                [item for sublist in new_SubId for item in sublist]
                            )
                        )  # flatten list and make entries unique

                    old_SubId = np.array(
                        [item for sublist in old_SubId for item in sublist],
                        dtype=np.int,
                    )  # flatten list
                    old_SubIds += list(old_SubId)

                old_SubIds = list(np.sort(np.unique(old_SubIds)))

                idx_basins = [
                    list(np.where(shape["SubId"] == oo)[0]) for oo in old_SubIds
                ]
                idx_basins = [
                    item for sublist in idx_basins for item in sublist
                ]  # flatten list
                idx_basins = list(
                    np.unique(idx_basins)
                )  # getting only unique list indexes

        else:  # all HRUs to be processed

            idx_basins = list(np.arange(0, len(shape)))

        # make sure HRUs are only once in this list
        hrus = np.array(shape.loc[idx_basins][key_colname])  # [sort_idx]

        idx_basins_unique = []
        hrus_unique = []
        for ihru, hru in enumerate(hrus):

            if not (hru in hrus_unique):

                hrus_unique.append(hrus[ihru])
                idx_basins_unique.append(idx_basins[ihru])

        idx_basins = idx_basins_unique
        hrus = hrus_unique

        # order according to values in "key_colname"; just to make sure outputs will be sorted in the end
        sort_idx = np.argsort(shape.loc[idx_basins][key_colname])
        print(
            "   >>> HRU_IDs found = ",
            list(np.array(shape.loc[idx_basins][key_colname])[sort_idx]),
            "  (total: ",
            len(idx_basins),
            ")",
        )

        # reduce the routing product dataset now to only what we will need
        shape = shape.loc[np.array(idx_basins)[sort_idx]]

        # indexes of all lines in df
        keys = shape.index
        nsubbasins = len(keys)

        # initialize
        coord_catch_wkt = {}

        # loop over all subbasins and transform coordinates into equal-area projection
        for kk in keys:

            ibasin = shape.loc[kk]

            poly = ibasin.geometry
            coord_catch_wkt[kk] = ogr.CreateGeometryFromWkt(poly.to_wkt())

        # -------------------------------
        # construct all grid cell polygons
        # -------------------------------
        print(" ")
        print("   (3) Generate shapes for NetCDF grid cells ...")

        grid_cell_geom_gpd_wkt = [[[] for ilon in range(nlon)] for ilat in range(nlat)]
        for ilat in range(nlat):
            if ilat % 10 == 0:
                print("   >>> Latitudes done: {0} of {1}".format(ilat, nlat))

            for ilon in range(nlon):

                # -------------------------
                # EPSG:3035   needs a swap before and after transform ...
                # -------------------------
                # gridcell_edges = [ [lath[ilat,ilon]    , lonh[ilat,  ilon]    ],            # for some reason need to switch lat/lon that transform works
                #                    [lath[ilat+1,ilon]  , lonh[ilat+1,ilon]    ],
                #                    [lath[ilat+1,ilon+1], lonh[ilat+1,ilon+1]  ],
                #                    [lath[ilat,ilon+1]  , lonh[ilat,  ilon+1]  ]]

                # tmp = shape_to_geometry(gridcell_edges, epsg=crs_caea)
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

                tmp = shape_to_geometry(gridcell_edges, epsg=crs_caea)
                grid_cell_geom_gpd_wkt[ilat][ilon] = tmp

        # -------------------------------
        # Derive overlay and calculate weights
        # -------------------------------
        print(" ")
        print("   (4) Deriving weights ...")

        filename = output_file
        ff = open(filename, "w")
        ff.write(":GridWeights                     \n")
        ff.write("   #                                \n")
        ff.write("   # [# HRUs]                       \n")
        ff.write("   :NumberHRUs       {0}            \n".format(nsubbasins))
        ff.write("   :NumberGridCells  {0}            \n".format(nlon * nlat))
        ff.write("   #                                \n")
        ff.write("   # [HRU ID] [Cell #] [w_kl]       \n")

        for ikk, kk in enumerate(keys):

            ibasin = shape.loc[kk]

            area_basin = coord_catch_wkt[kk].Area()
            enve_basin = coord_catch_wkt[
                kk
            ].GetEnvelope()  # bounding box around basin (for easy check of proximity)

            area_all = 0.0
            ncells = 0

            for ilat in range(nlat):
                for ilon in range(nlon):

                    enve_gridcell = grid_cell_geom_gpd_wkt[ilat][
                        ilon
                    ].GetEnvelope()  # bounding box around grid-cell (for easy check of proximity)
                    grid_is_close = check_proximity_of_envelops(
                        enve_gridcell, enve_basin
                    )

                    if (
                        grid_is_close
                    ):  # this check decreases runtime DRASTICALLY (from ~6h to ~1min)

                        grid_cell_area = grid_cell_geom_gpd_wkt[ilat][ilon].Area()

                        inter = grid_cell_geom_gpd_wkt[ilat][ilon].Intersection(
                            coord_catch_wkt[kk].Buffer(0.0)
                        )  # "fake" buffer to avoid invalid polygons and weirdos dumped by ArcGIS
                        area_intersect = inter.Area()

                        area_all += area_intersect
                        if area_intersect > 0:
                            ncells += 1

                            print(
                                "   >>> {0},{1},{2},{3},{4}".format(
                                    int(ibasin[key_colname]),
                                    ilat,
                                    ilon,
                                    ilat * nlon + ilon,
                                    area_intersect / area_basin,
                                )
                            )
                            ff.write(
                                "   {0}   {1}   {2}\n".format(
                                    int(ibasin[key_colname]),
                                    ilat * nlon + ilon,
                                    area_intersect / area_basin,
                                )
                            )

            print(
                "   >>> (Sub-)Basin: {0} ({1} of {2})".format(
                    int(ibasin[key_colname]), ikk + 1, nsubbasins
                )
            )
            print("   >>> Derived area of {0}  cells: {1}".format(ncells, area_all))
            print("   >>> Read area from shapefile:   {0}".format(area_basin))
            print(
                "   >>> error:                      {0}%".format(
                    (area_basin - area_all) / area_basin * 100.0
                )
            )
            print("   ")

        ff.write(":EndGridWeights \n")
        ff.close()

        print("")
        print("Wrote: ", filename)
        print("")

    def create_gridcells_from_centers(self, lat, lon):

        # create array of edges where (x,y) are always center cells
        nlon = np.shape(lon)[1]
        nlat = np.shape(lat)[0]
        lonh = np.empty((nlat + 1, nlon + 1), dtype=np.float)
        lath = np.empty((nlat + 1, nlon + 1), dtype=np.float)
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

    def shape_to_geometry(self, shape_from_jsonfile, epsg=None):

        # converts shape read from shapefile to geometry
        # epsg :: integer EPSG code

        ring_shape = ogr.Geometry(ogr.wkbLinearRing)

        for ii in shape_from_jsonfile:
            ring_shape.AddPoint_2D(ii[0], ii[1])
        # close ring
        ring_shape.AddPoint_2D(shape_from_jsonfile[0][0], shape_from_jsonfile[0][1])

        poly_shape = ogr.Geometry(ogr.wkbPolygon)
        poly_shape.AddGeometry(ring_shape)

        if not (epsg is None):
            source = osr.SpatialReference()
            source.ImportFromEPSG(crs_lldeg)  # usual lat/lon projection

            target = osr.SpatialReference()
            target.ImportFromEPSG(epsg)  # any projection to convert to

            transform = osr.CoordinateTransformation(source, target)
            poly_shape.Transform(transform)

        return poly_shape

    def check_proximity_of_envelops(self, gridcell_envelop, shape_envelop):

        # checks if two envelops are in proximity (intersect)

        # minX  --> env[0]
        # maxX  --> env[1]
        # minY  --> env[2]
        # maxY  --> env[3]

        if (
            (gridcell_envelop[0] <= shape_envelop[1])
            and (gridcell_envelop[1] >= shape_envelop[0])
            and (gridcell_envelop[2] <= shape_envelop[3])
            and (gridcell_envelop[3] >= shape_envelop[2])
        ):

            grid_is_close = True

        else:

            grid_is_close = False

        return grid_is_close

    def check_gridcell_in_proximity_of_shape(self, gridcell_edges, shape_from_jsonfile):

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

        if (
            (min_lat_cell <= max_lat_shape)
            and (max_lat_cell >= min_lat_shape)
            and (min_lon_cell <= max_lon_shape)
            and (max_lon_cell >= min_lon_shape)
        ):

            grid_is_close = True

        else:

            grid_is_close = False

        return grid_is_close
