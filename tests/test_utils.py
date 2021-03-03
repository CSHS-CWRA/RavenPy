import tempfile
from pathlib import Path

import numpy as np
import pytest


from ravenpy.utilities.testdata import get_local_testdata


class TestOperations:
    utils = pytest.importorskip("ravenpy.utils")

    zipped_file = get_local_testdata("polygons/mars.zip")
    non_zipped_file = get_local_testdata("polygons/mars.geojson")

    def test_circular_mean_aspect(self):
        northern_angles = np.array([330, 30, 15, 345])
        slight_northeast_angles = np.append(northern_angles, [0.000001])
        eastern_angles = np.arange(45, 125, 1.25)
        southwest_angles = np.array([181, 182.25, 183.5, 222])

        assert self.utils.circular_mean_aspect(northern_angles) == 360
        np.testing.assert_almost_equal(
            self.utils.circular_mean_aspect(slight_northeast_angles), 0, decimal=3
        )
        assert self.utils.circular_mean_aspect(eastern_angles) == 84.375
        np.testing.assert_almost_equal(
            self.utils.circular_mean_aspect(southwest_angles), 191.88055987
        )

    def test_parse_lonlat(self):
        lonlats = ["123,345", "122233344.11111 34554554.2", "1111,2222"]
        for ll in lonlats:
            lonlat = self.utils.parse_lonlat(ll)
            assert len(lonlat) == 2
            assert isinstance(lonlat[0], float)
            assert isinstance(lonlat[1], float)

        with pytest.raises(Exception):
            self.utils.parse_lonlat("This isn't a number, 333.444")

    def test_address_append(self):
        non_existing_tarred_file = "polygons.tar"

        assert "zip://" in self.utils.address_append(self.zipped_file)
        assert "tar://" in self.utils.address_append(non_existing_tarred_file)
        assert not self.utils.address_append(self.non_zipped_file).startswith(("zip://", "tar://"))

    def test_archive_sniffer(self):
        probable_shp = self.utils.archive_sniffer(self.zipped_file)
        assert probable_shp == ["/tmp/mars.shp"]

        probable_shp = self.utils.archive_sniffer(self.zipped_file, working_dir="/tmp")
        assert probable_shp == ["/tmp/mars.shp"]

    def test_archive_extract(self):

        assert self.zipped_file.exists()

        files = list()
        with tempfile.TemporaryDirectory() as tdir:
            files.extend(self.utils.generic_extract_archive(self.zipped_file, output_dir=tdir))
            assert len(files) == 5
            for f in files:
                assert Path(f).exists()
        assert not np.any([Path(f).exists() for f in files])

        files = self.utils.generic_extract_archive(self.zipped_file)
        assert np.all([Path(f).exists() for f in files])


class TestFileInfoFuncs:
    utils = pytest.importorskip("ravenpy.utils")

    zipped_file = get_local_testdata("polygons/mars.zip")
    geojson_file = get_local_testdata("polygons/mars.geojson")
    raster_file = get_local_testdata("nasa/Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff")

    non_existing_file = "unreal.zip"

    def test_raster_datatype_sniffer(self):
        datatype = self.utils.raster_datatype_sniffer(self.raster_file)
        assert datatype.lower() == "uint8"

    def test_crs_sniffer(self):
        assert self.utils.crs_sniffer(self.zipped_file) == 4326
        assert set(self.utils.crs_sniffer(self.geojson_file, self.raster_file)) == {4326}

    def test_single_file_check(self):
        one = [Path(__file__).parent / "__init__.py"]
        zero = []
        three = [1, Path().root, 2.333]

        assert self.utils.single_file_check(one) == one[0]

        with pytest.raises(FileNotFoundError):
            self.utils.single_file_check(zero)

        with pytest.raises(NotImplementedError):
            self.utils.single_file_check(three)

    def test_boundary_check(self):
        # NOTE: does not presently accept zipped files.

        with pytest.warns(None):
            self.utils.boundary_check([self.geojson_file, self.raster_file], max_y=80)

        with pytest.warns(UserWarning):
            self.utils.boundary_check([self.geojson_file, self.raster_file], max_y=15)

        with pytest.raises(FileNotFoundError):
            self.utils.boundary_check([self.non_existing_file])

    @pytest.mark.skip(reason="Not presently testable")
    def test_multipolygon_check(self):
        pass


class TestGdalOgrFunctions:
    utils = pytest.importorskip("ravenpy.utils")
    fiona = pytest.importorskip("fiona")
    sgeo = pytest.importorskip("shapely.geometry")

    geojson_file = get_local_testdata("polygons/mars.geojson")
    raster_file = get_local_testdata("nasa/Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff")

    def test_gdal_aspect_not_projected(self):
        aspect_grid = self.utils.gdal_aspect_analysis(self.raster_file)
        np.testing.assert_almost_equal(
            self.utils.circular_mean_aspect(aspect_grid), 10.9119033
        )

        # test with creation of a temporary file
        aspect_tempfile = tempfile.NamedTemporaryFile(
            prefix="aspect_", suffix=".tiff", delete=False
        ).name
        aspect_grid = self.utils.gdal_aspect_analysis(
            self.raster_file, set_output=aspect_tempfile
        )
        np.testing.assert_almost_equal(
            self.utils.circular_mean_aspect(aspect_grid), 10.9119033
        )
        assert Path(aspect_tempfile).stat().st_size > 0

    # Slope values are high due to data values using Geographic CRS
    def test_gdal_slope_not_projected(self):
        slope_grid = self.utils.gdal_slope_analysis(self.raster_file)
        np.testing.assert_almost_equal(slope_grid.min(), 0.0)
        np.testing.assert_almost_equal(slope_grid.mean(), 64.4365427)
        np.testing.assert_almost_equal(slope_grid.max(), 89.71747, 5)

        slope_tempfile = tempfile.NamedTemporaryFile(
            prefix="slope_", suffix=".tiff", delete=False
        ).name
        slope_grid = self.utils.gdal_slope_analysis(self.raster_file, set_output=slope_tempfile)
        np.testing.assert_almost_equal(slope_grid.mean(), 64.4365427)
        assert Path(slope_tempfile).stat().st_size > 0

    # Slope values are high due to data values using Geographic CRS
    def test_dem_properties(self):
        dem_properties = self.utils.dem_prop(self.raster_file)
        np.testing.assert_almost_equal(dem_properties["aspect"], 10.9119033)
        np.testing.assert_almost_equal(dem_properties["elevation"], 79.0341721)
        np.testing.assert_almost_equal(dem_properties["slope"], 64.4365427)

        with self.fiona.open(self.geojson_file) as gj:
            feature = next(iter(gj))
            geom = self.sgeo.shape(feature["geometry"])

        region_dem_properties = self.utils.dem_prop(self.raster_file, geom=geom)
        np.testing.assert_almost_equal(region_dem_properties["aspect"], 280.6814208)
        np.testing.assert_almost_equal(region_dem_properties["elevation"], 145.8899082)
        np.testing.assert_almost_equal(region_dem_properties["slope"], 61.2650882)

    # Slope values are high due to data values using Geographic CRS
    def test_geom_properties(self):
        with self.fiona.open(self.geojson_file) as gj:
            iterable = iter(gj)
            feature_1 = next(iterable)
            feature_2 = next(iterable)
            geom_1 = self.sgeo.shape(feature_1["geometry"])
            geom_2 = self.sgeo.shape(feature_2["geometry"])

        geom_1_properties = self.utils.geom_prop(geom_1)
        np.testing.assert_almost_equal(geom_1_properties["area"], 357.9811899)
        np.testing.assert_almost_equal(
            geom_1_properties["centroid"], (-128.3959836, 19.1572278)
        )
        np.testing.assert_almost_equal(geom_1_properties["perimeter"], 68.4580077)
        np.testing.assert_almost_equal(geom_1_properties["gravelius"], 1.0206790)

        geom_2_properties = self.utils.geom_prop(geom_2)
        np.testing.assert_almost_equal(geom_2_properties["area"], 361.5114221)
        np.testing.assert_almost_equal(
            geom_2_properties["centroid"], (-70.2394629, 45.7698029)
        )
        np.testing.assert_almost_equal(geom_2_properties["perimeter"], 96.1035859)
        np.testing.assert_almost_equal(geom_2_properties["gravelius"], 1.4258493)


class TestGenericGeoOperations:
    utils = pytest.importorskip("ravenpy.utils")
    fiona = pytest.importorskip("fiona")
    rasterio = pytest.importorskip("rasterio")
    sgeo = pytest.importorskip("shapely.geometry")

    geojson_file = get_local_testdata("polygons/mars.geojson")
    raster_file = get_local_testdata("nasa/Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff")

    def test_vector_reprojection(self):
        # TODO: It would be awesome if this returned a temporary filepath if no file given.
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".geojson", delete=False
        ).name
        self.utils.generic_vector_reproject(
            self.geojson_file, projected=reproj_file, target_crs="EPSG:3348"
        )

        with self.fiona.open(reproj_file) as gj:
            iterable = iter(gj)
            feature = next(iterable)
            geom = self.sgeo.shape(feature["geometry"])

        geom_properties = self.utils.geom_prop(geom)
        np.testing.assert_almost_equal(geom_properties["area"], 6450001762792.884, 3)
        np.testing.assert_almost_equal(
            geom_properties["centroid"], (1645777.7589835, -933242.1203143)
        )
        np.testing.assert_almost_equal(geom_properties["perimeter"], 9194343.1759303)
        np.testing.assert_almost_equal(geom_properties["gravelius"], 1.0212589)

    def test_raster_warp(self):
        # TODO: It would be awesome if this returned a temporary filepath if no file given.
        # TODO: either use `output` or `reprojected/warped` for these functions.
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".tiff", delete=False
        ).name
        self.utils.generic_raster_warp(
            self.raster_file, output=reproj_file, target_crs="EPSG:3348"
        )

        with self.rasterio.open(reproj_file) as gt:
            assert gt.crs.to_epsg() == 3348
            np.testing.assert_almost_equal(gt.bounds.left, -2077535.25979486)
            np.testing.assert_almost_equal(gt.bounds.right, 15591620.75098695)
            np.testing.assert_almost_equal(gt.bounds.bottom, -4167898.76317739)
            np.testing.assert_almost_equal(gt.bounds.top, 5817014.91999878)

            data = gt.read(1)  # read band 1 (red)
            assert data.min() == 0
            assert data.max() == 255
            np.testing.assert_almost_equal(data.mean(), 60.7291936)

    def test_warped_raster_slope(self):
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".tiff", delete=False
        ).name
        self.utils.generic_raster_warp(
            self.raster_file, output=reproj_file, target_crs="EPSG:3348"
        )
        slope_grid = self.utils.gdal_slope_analysis(reproj_file)

        np.testing.assert_almost_equal(slope_grid.min(), 0.0)
        np.testing.assert_almost_equal(slope_grid.mean(), 0.0034991)
        np.testing.assert_almost_equal(slope_grid.max(), 0.3523546)

    def test_warped_raster_aspect(self):
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".tiff", delete=False
        ).name
        self.utils.generic_raster_warp(
            self.raster_file, output=reproj_file, target_crs="EPSG:3348"
        )
        aspect_grid = self.utils.gdal_aspect_analysis(reproj_file)

        np.testing.assert_almost_equal(
            self.utils.circular_mean_aspect(aspect_grid), 7.7805879
        )

    def test_raster_clip(self):
        with self.fiona.open(self.geojson_file) as gj:
            feature = next(iter(gj))
            geom = self.sgeo.shape(feature["geometry"])

        clipped_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".tiff", delete=False
        ).name
        self.utils.generic_raster_clip(self.raster_file, clipped_file, geometry=geom)

        with self.rasterio.open(clipped_file) as gt:
            assert gt.crs.to_epsg() == 4326

            data = gt.read(1)  # read band 1 (red)
            assert data.min() == 0
            assert data.max() == 255
            np.testing.assert_almost_equal(data.mean(), 102.8222965)

    def test_shapely_pyproj_transform(self):
        with self.fiona.open(self.geojson_file) as gj:
            feature = next(iter(gj))
            geom = self.sgeo.shape(feature["geometry"])

        transformed = self.utils.geom_transform(geom, target_crs="EPSG:3348")
        np.testing.assert_almost_equal(
            transformed.bounds,
            (188140.3820599, -2374936.1363096, 3086554.0207066, 409691.2180337),
        )
        np.testing.assert_almost_equal(transformed.centroid.x, 1645777.7589835)
        np.testing.assert_almost_equal(transformed.centroid.y, -933242.1203143)
        np.testing.assert_almost_equal(transformed.area, 6450001762792.884, 3)


class TestGIS:
    utils = pytest.importorskip("ravenpy.utils")
    sgeo = pytest.importorskip("shapely.geometry")

    vector_file = get_local_testdata("polygons/mars.geojson")

    def test_get_bbox_single(self):
        w, s, n, e = self.utils.get_bbox(self.vector_file, all_features=False)
        np.testing.assert_almost_equal(w, -139.8514262)
        np.testing.assert_almost_equal(s, 8.3754794)
        np.testing.assert_almost_equal(n, -117.4753973)
        np.testing.assert_almost_equal(e, 29.6327068)

    def test_get_bbox_all(self):
        w, s, n, e = self.utils.get_bbox(self.vector_file)
        np.testing.assert_almost_equal(w, -139.8514262)
        np.testing.assert_almost_equal(s, 8.3754794)
        np.testing.assert_almost_equal(n, -38.7397456)
        np.testing.assert_almost_equal(e, 64.1757015)

    def test_feature_contains(self):
        point = -69.0, 45
        assert isinstance(self.utils.feature_contains(point, self.vector_file), dict)
        assert isinstance(self.utils.feature_contains(self.sgeo.Point(point), self.vector_file), dict)
