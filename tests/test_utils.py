import pytest

from ravenpy import RavenPyDependencyError

try:
    import ravenpy.utils as utils
except RavenPyDependencyError:
    utils = False


class TestOperations:

    def test_circular_mean_aspect(self):
        pass

    def test_parse_lonlat(self):
        pass

    def test_address_append(self):
        pass

    def test_archive_extract(self):
        pass


@pytest.mark.skipif(condition=utils is False, reason="GIS dependencies are needed.")
class TestFileInfoFuncs:

    def test_raster_sniffers(self):
        pass

    def test_vector_sniffers(self):
        pass

    def test_archive_sniffer(self):
        pass

    def test_crs_sniffer(self):
        pass

    def test_invalid_crs_sniffer(self):
        pass

    def test_single_file_check(self):
        pass

    def test_boundary_check(self):
        pass

    def test_multipolygon_check(self):
        pass


@pytest.mark.skipif(condition=utils is False, reason="GIS dependencies are needed.")
class TestGdalOgrFunctions:

    def test_gdal_aspect(self):
        pass

    def test_gdal_slope(self):
        pass

    def test_dem_properties(self):
        pass

    def test_geom_properties(self):
        pass


@pytest.mark.skipif(condition=utils is False, reason="GIS dependencies are needed.")
class TestGeoOperations:

    def test_vector_reprojection(self):
        pass

    def test_raster_warp(self):
        pass

    def test_raster_clip(self):
        pass

    def test_shapely_pyproj_transform(self):
        pass
