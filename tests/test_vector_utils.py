import tempfile

import geojson
import numpy as np
import pytest
import shapely.geometry as sgeo

from ravenpy.utilities import vector
from ravenpy.utilities.testdata import get_local_testdata


class TestVectorUtils:
    geojson_file = get_local_testdata("polygons/mars.geojson")
    routing_product_shapefile = get_local_testdata(
        "raven-routing-sample/finalcat_hru_info.zip"
    )

    def test_shapely_pyproj_geojson_transform(self):
        feat = geojson.load(open(self.geojson_file))
        geom = sgeo.shape(feat[0].geometry)

        transformed = vector.geom_transform(geom, target_crs="EPSG:3348")
        np.testing.assert_almost_equal(
            transformed.bounds, (188140, -2374936, 3086554, 409691), 0
        )
        np.testing.assert_almost_equal(transformed.centroid.x, 1645777, 0)
        np.testing.assert_almost_equal(transformed.centroid.y, -933242, 0)
        np.testing.assert_almost_equal(transformed.area, 6450001868342, 0)

    @pytest.mark.slow
    def test_shapely_pyproj_shapely_transform_properties(self, tmp_path):
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".geojson", delete=False, dir=tmp_path
        ).name
        vector.generic_vector_reproject(
            self.routing_product_shapefile,
            projected=reproj_file,
            target_crs="EPSG:3348",
        )

        feat = geojson.load(open(reproj_file))
        geom = sgeo.shape(feat[0].geometry)

        geom_properties = vector.geom_prop(geom)
        np.testing.assert_almost_equal(geom_properties["area"], 310646674, 0)
        np.testing.assert_almost_equal(
            geom_properties["centroid"], (7460914.0, 1413210.6), 1
        )
        np.testing.assert_almost_equal(geom_properties["perimeter"], 371435.385, 3)
        np.testing.assert_almost_equal(geom_properties["gravelius"], 5.9449059)

    def test_shapely_pyproj_geojson_transform_properties(self, tmp_path):
        # TODO: It would be awesome if this returned a temporary filepath if no file given.
        reproj_file = tempfile.NamedTemporaryFile(
            prefix="reproj_", suffix=".geojson", delete=False, dir=tmp_path
        ).name
        vector.generic_vector_reproject(
            self.geojson_file, projected=reproj_file, target_crs="EPSG:3348"
        )

        feat = geojson.load(open(reproj_file))
        geom = sgeo.shape(feat[0].geometry)

        geom_properties = vector.geom_prop(geom)
        np.testing.assert_almost_equal(geom_properties["area"], 6450001868342, 0)
        np.testing.assert_almost_equal(
            geom_properties["centroid"], (1645777.7, -933242.1), 1
        )
        np.testing.assert_almost_equal(geom_properties["perimeter"], 9194343, 0)
        np.testing.assert_almost_equal(geom_properties["gravelius"], 1.0212589)

    # Slope values are high due to data values using Geographic CRS
    def test_geom_properties_multi_feature(self):
        feat = geojson.load(open(self.geojson_file))
        geom_0 = sgeo.shape(feat[0].geometry)
        geom_1 = sgeo.shape(feat[1].geometry)

        geom_1_properties = vector.geom_prop(geom_0)
        np.testing.assert_almost_equal(geom_1_properties["area"], 357.981, 3)
        np.testing.assert_almost_equal(
            geom_1_properties["centroid"], (-128.396, 19.157), 3
        )
        np.testing.assert_almost_equal(geom_1_properties["perimeter"], 68.458, 3)
        np.testing.assert_almost_equal(geom_1_properties["gravelius"], 1.021, 3)

        geom_2_properties = vector.geom_prop(geom_1)
        np.testing.assert_almost_equal(geom_2_properties["area"], 361.511, 3)
        np.testing.assert_almost_equal(
            geom_2_properties["centroid"], (-70.239, 45.770), 3
        )
        np.testing.assert_almost_equal(geom_2_properties["perimeter"], 96.104, 3)
        np.testing.assert_almost_equal(geom_2_properties["gravelius"], 1.426, 3)
