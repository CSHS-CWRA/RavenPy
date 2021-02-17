import pytest
import tempfile
from pathlib import Path

import numpy as np

from ravenpy import RavenPyDependencyError

try:
    import ravenpy.utils as utils
    import fiona
    import rasterio
    from shapely.geometry import shape, GeometryCollection
except RavenPyDependencyError:
    utils = False

from .common import test_data


class TestOperations:
    def test_circular_mean_aspect(self):
        northern_angles = np.array([330, 30, 15, 345])
        slight_northeast_angles = np.append(northern_angles, [0.000001])
        eastern_angles = np.arange(45, 125, 1.25)
        southwest_angles = np.array([181, 182.25, 183.5, 222])

        assert utils.circular_mean_aspect(northern_angles) == 360
        np.testing.assert_almost_equal(
            utils.circular_mean_aspect(slight_northeast_angles), 0, decimal=3
        )
        assert utils.circular_mean_aspect(eastern_angles) == 84.375
        np.testing.assert_almost_equal(
            utils.circular_mean_aspect(southwest_angles), 191.88055987
        )

    def test_parse_lonlat(self):
        lonlats = ["123,345", "122233344.11111 34554554.2", "1111,2222"]
        for ll in lonlats:
            lonlat = utils.parse_lonlat(ll)
            assert len(lonlat) == 2
            assert isinstance(lonlat[0], float)
            assert isinstance(lonlat[1], float)

        with pytest.raises(Exception):
            utils.parse_lonlat("This isn't a number, 333.444")

    def test_address_append(self):
        zipped_file = test_data() / "polygons.zip"
        tarred_file = test_data() / "polygons.tar"
        non_zipped_file = test_data() / "polygons.geojson"

        assert "zip://" in utils.address_append(zipped_file)
        assert "tar://" in utils.address_append(tarred_file)
        # Need to change return type in RAVEN address_append to always be str
        assert not str(utils.address_append(non_zipped_file)).startswith(
            ("zip://", "tar://")
        )

    def test_archive_sniffer(self):
        zipped_file = test_data() / "polygons.zip"

        # `working_dir` should be allowed as None
        probable_shp = utils.archive_sniffer(zipped_file, working_dir="/tmp")
        assert probable_shp == ["/tmp/polygons.shp"]

    def test_archive_extract(self):
        zipped_file = test_data() / "polygons.zip"
        assert zipped_file.exists()

        files = list()
        with tempfile.TemporaryDirectory() as tdir:
            files.extend(utils.generic_extract_archive(zipped_file, output_dir=tdir))
            assert len(files) == 5
            for f in files:
                assert Path(f).exists()
        assert not np.any([Path(f).exists() for f in files])

        files = utils.generic_extract_archive(zipped_file)
        assert np.all([Path(f).exists() for f in files])


@pytest.mark.skipif(condition=utils is False, reason="GIS dependencies are needed.")
class TestFileInfoFuncs:
    def test_raster_datatype_sniffer(self):
        raster_file = (
            test_data() / "Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff"
        )

        datatype = utils.raster_datatype_sniffer(raster_file)
        assert datatype.lower() == "uint8"

    def test_crs_sniffer(self):
        # FIXME: This utility should not complain if given a single file / list
        zipped_file = test_data() / "polygons.zip"
        geojson_file = test_data() / "polygons.geojson"
        raster_file = (
            test_data() / "Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff"
        )

        # FIXME: Should this be raising a FileNotFound internally? Probably not.
        with pytest.raises(Exception):
            utils.crs_sniffer(zipped_file)

        # TODO: This will fail with the new PyProj when ported. Will be == int(4326).
        assert set(utils.crs_sniffer(geojson_file, raster_file)) == {"+init=epsg:4326"}

    @pytest.mark.skip
    def test_single_file_check(self):
        # FIXME: This utility should ensure that files exist. Everything goes right now.
        # FIXME: Exceptions do not work as intended. Port changes from RAVEN.
        one = [Path(__file__).parent / "__init__.py"]
        zero = []
        three = [1, Path().root, 2.333]

        assert utils.single_file_check(one) == one[0]

        with pytest.raises(FileNotFoundError):
            utils.single_file_check(zero)

        with pytest.raises(NotImplementedError):
            utils.single_file_check(three)

    @pytest.mark.skip
    def test_boundary_check(self):
        # FIXME: This utility should not complain if given a single file / list
        # FIXME: This is very broken. Needs to be fully rewritten.
        zipped_file = test_data() / "polygons.zip"
        geojson_file = test_data() / "polygons.geojson"
        raster_file = (
            test_data() / "Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff"
        )
        nonexists_file = test_data() / "unreal.zip"

        with pytest.warns(None):
            utils.boundary_check([zipped_file, geojson_file, raster_file], max_y=80)

        with pytest.warns(UserWarning):
            utils.boundary_check([zipped_file, geojson_file, raster_file], max_y=15)

        with pytest.raises(FileNotFoundError):
            utils.boundary_check([nonexists_file])

    @pytest.mark.skip(reason="Not presently testable")
    def test_multipolygon_check(self):
        pass


@pytest.mark.skipif(condition=utils is False, reason="GIS dependencies are needed.")
class TestGdalOgrFunctions:

    raster_file = test_data() / "Mars_MGS_MOLA_DEM_georeferenced_region_compressed.tiff"
    geojson_file = test_data() / "polygons.geojson"

    def test_gdal_aspect_not_projected(self):
        # FIXME: This should remove the temporary file, saving the grid in memory only.
        aspect_grid = utils.gdal_aspect_analysis(self.raster_file)
        np.testing.assert_almost_equal(
            utils.circular_mean_aspect(aspect_grid), 10.9119033
        )

        # test with creation of a temporary file
        aspect_tempfile = tempfile.NamedTemporaryFile(
            prefix="aspect_", suffix=".tiff", delete=False
        ).name
        aspect_grid = utils.gdal_aspect_analysis(
            self.raster_file, output=aspect_tempfile
        )
        np.testing.assert_almost_equal(
            utils.circular_mean_aspect(aspect_grid), 10.9119033
        )
        assert Path(aspect_tempfile).stat().st_size > 0

    # Slope values are high due to data values using Geographic CRS
    def test_gdal_slope_not_projected(self):
        # FIXME: This should remove the temporary file, saving the grid in memory only.
        slope_grid = utils.gdal_slope_analysis(self.raster_file)
        np.testing.assert_almost_equal(slope_grid.min(), 0.0)
        np.testing.assert_almost_equal(slope_grid.mean(), 64.4365427)
        np.testing.assert_almost_equal(slope_grid.max(), 89.71747, 5)

        slope_tempfile = tempfile.NamedTemporaryFile(
            prefix="slope_", suffix=".tiff", delete=False
        ).name
        slope_grid = utils.gdal_slope_analysis(self.raster_file, output=slope_tempfile)
        np.testing.assert_almost_equal(slope_grid.mean(), 64.4365427)
        assert Path(slope_tempfile).stat().st_size > 0

    # Slope values are high due to data values using Geographic CRS
    def test_dem_properties(self):
        dem_properties = utils.dem_prop(self.raster_file)
        np.testing.assert_almost_equal(dem_properties["aspect"], 10.9119033)
        np.testing.assert_almost_equal(dem_properties["elevation"], 79.0341721)
        np.testing.assert_almost_equal(dem_properties["slope"], 64.4365427)

        with fiona.open(self.geojson_file) as gj:
            feature = next(iter(gj))
            geom = shape(feature["geometry"])

        region_dem_properties = utils.dem_prop(self.raster_file, geom=geom)
        np.testing.assert_almost_equal(region_dem_properties["aspect"], 280.6814208)
        np.testing.assert_almost_equal(region_dem_properties["elevation"], 145.8899082)
        np.testing.assert_almost_equal(region_dem_properties["slope"], 61.2650882)

    # Slope values are high due to data values using Geographic CRS
    def test_geom_properties(self):
        with fiona.open(self.geojson_file) as gj:
            iterable = iter(gj)
            feature_1 = next(iterable)
            feature_2 = next(iterable)
            geom_1 = shape(feature_1["geometry"])
            geom_2 = shape(feature_2["geometry"])

        geom_1_properties = utils.geom_prop(geom_1)
        np.testing.assert_almost_equal(geom_1_properties["area"], 357.9811899)
        np.testing.assert_almost_equal(
            geom_1_properties["centroid"], (-128.3959836, 19.1572278)
        )
        np.testing.assert_almost_equal(geom_1_properties["perimeter"], 68.4580077)
        np.testing.assert_almost_equal(geom_1_properties["gravelius"], 1.0206790)

        geom_2_properties = utils.geom_prop(geom_2)
        np.testing.assert_almost_equal(geom_2_properties["area"], 361.5114221)
        np.testing.assert_almost_equal(
            geom_2_properties["centroid"], (-70.2394629,  45.7698029)
        )
        np.testing.assert_almost_equal(geom_2_properties["perimeter"], 96.1035859)
        np.testing.assert_almost_equal(geom_2_properties["gravelius"], 1.4258493)


@pytest.mark.skipif(condition=utils is False, reason="GIS dependencies are needed.")
class TestGeoOperations:
    def test_vector_reprojection(self):
        pass

    def test_raster_warp(self):
        pass

    def test_warped_raster_slope(self):
        pass

    def test_warped_raster_aspect(self):
        pass

    def test_raster_clip(self):
        pass

    def test_shapely_pyproj_transform(self):
        pass
