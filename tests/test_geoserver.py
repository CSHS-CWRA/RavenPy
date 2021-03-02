import tempfile

import numpy as np
import pytest
import ravenpy.utils

try:
    import fiona
    import geopandas as gpd
    import rasterio
    from shapely.geometry import GeometryCollection, shape, Point

    from ravenpy.utilities import geoserver
    import ravenpy.utils as utils
except (ModuleNotFoundError, ImportError):
    utils = False

from ravenpy.utilities.testdata import get_local_testdata


class TestHydroBASINS:
    def test_select_hybas_na_domain(self):
        bbox = (-68.0, 50.0) * 2
        dom = geoserver.select_hybas_domain(bbox)
        assert dom == "na"

    def test_select_hybas_ar_domain(self):
        bbox = (-114.65, 61.35) * 2
        dom = geoserver.select_hybas_domain(bbox)
        assert dom == "ar"

    def test_get_hydrobasins_location_wfs(self):
        lake_winnipeg = (-98.03575958286369, 52.88238524279493)
        feature = geoserver.get_hydrobasins_location_wfs(
            coordinates=lake_winnipeg * 2, lakes=True, domain="na"
        )

        gml = tempfile.NamedTemporaryFile(suffix=".gml", delete=False)
        with open(gml.name, "wb") as f:
            f.write(feature)

        with fiona.open(gml.name) as src:
            iterable = iter(src)
            feat = next(iterable)
            geom = shape(feat["geometry"])
            assert geom.bounds == (-99.2731, 50.3603, -96.2578, 53.8705)
            np.testing.assert_almost_equal(geom.area, 3.2530867)

    def test_get_hydrobasins_attributes_wfs(self):
        rio_grande = (-80.475, 8.4)
        feature = geoserver.get_hydrobasins_location_wfs(
            coordinates=rio_grande * 2, lakes=True, domain="na"
        )

        gml = tempfile.NamedTemporaryFile(suffix=".gml", delete=False)
        with open(gml.name, "wb") as f:
            f.write(feature)

        with fiona.open(gml.name) as src:
            iterable = iter(src)
            feat = next(iterable)
            main_bas = feat["properties"]["MAIN_BAS"]

        # TODO: It would be swell to just have this function determine the domain if not given.
        region_url = geoserver.get_hydrobasins_attributes_wfs(
            attribute="MAIN_BAS", value=main_bas, domain="na"
        )
        gdf = gpd.read_file(region_url)

        assert len(gdf) == 18
        assert gdf.crs.to_epsg() == 4326
        assert set(gdf["MAIN_BAS"].values) == {7120000210}

        basin = gdf.dissolve(by="MAIN_BAS")
        np.testing.assert_equal(
            basin.geometry.bounds.values,
            np.array([[-80.8542, 8.2459, -80.1375, 8.7004]]),
        )

    def test_hydrobasins_upstream_ids_aggregate(self):
        puerto_cortes = (-83.525, 8.96)
        feature = geoserver.get_hydrobasins_location_wfs(
            coordinates=puerto_cortes * 2, lakes=True, domain="na"
        )

        gml = tempfile.NamedTemporaryFile(suffix=".gml", delete=False)
        with open(gml.name, "wb") as f:
            f.write(feature)

        with fiona.open(gml.name) as src:
            iterable = iter(src)
            feat = next(iterable)
            main_bas = feat["properties"]["MAIN_BAS"]
            hybas_id = feat["properties"]["HYBAS_ID"]

        # TODO: It would be swell to just have this function determine the domain if not given.
        region_url = geoserver.get_hydrobasins_attributes_wfs(
            attribute="MAIN_BAS", value=main_bas, domain="na"
        )
        gdf = gpd.read_file(region_url)
        gdf_upstream = geoserver.hydrobasins_upstream_ids(hybas_id, gdf)
        assert len(gdf) == len(gdf_upstream) + 1

        # FIXME: This file write step workaround is needed for some unknown reason.
        with tempfile.NamedTemporaryFile(prefix="hybas_", suffix=".json") as tf:
            gdf_upstream.to_file(tf.name, driver="GeoJSON")
            gdf_upstream = gpd.read_file(tf.name)
        aggregated = geoserver.hydrobasins_aggregate(gdf_upstream)

        assert len(aggregated) == 1
        assert aggregated.SUB_AREA.values == 4977.8
        np.testing.assert_equal(
            aggregated.geometry.bounds.values,
            np.array([[-83.8167, 8.7625, -82.7125, 9.5875]]),
        )


@pytest.mark.skip(reason="Not yet written.")
class TestHydroRouting:
    @pytest.mark.skip(reason="This function doesn't exist yet.")
    def test_hydro_routing_aggregate(self):
        pass

    @pytest.mark.skip(reason="This function doesn't exist yet.")
    def test_hydro_routing_upstream_ids(self):
        pass


@pytest.mark.skip(reason="Unsure if these tests are necessary.")
class TestWFS:
    def test_get_location_wfs(self):
        pass

    def test_get_feature_attributes_wfs(self):
        pass


class TestWCS:

    vector_file = get_local_testdata("polygons/Saskatoon.geojson")

    def test_get_raster_wcs(self):
        # TODO: This CRS needs to be redefined using modern pyproj-compatible strings.
        nalcms_crs = '+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs=True'

        with tempfile.NamedTemporaryFile(prefix="reprojected_", suffix=".json") as projected:
            utils.generic_vector_reproject(self.vector_file, projected.name, target_crs=nalcms_crs)
            bbox = ravenpy.utils.get_bbox(projected.name)

        raster_url = "public:CEC_NALCMS_LandUse_2010"
        raster_bytes = geoserver.get_raster_wcs(bbox, geographic=False, layer=raster_url)

        with tempfile.NamedTemporaryFile(prefix="wcs_", suffix=".tiff") as raster_file:
            with open(raster_file.name, "wb") as f:
                f.write(raster_bytes)

            with rasterio.open(raster_file) as src:
                assert src.width == 650
                assert src.height == 745
                np.testing.assert_array_equal(src.lnglat(), (-106.64193764047552, 52.1564202369763))

                data = src.read()
                assert np.unique(data).tolist() == [1, 5, 8, 10, 14, 15, 16, 17, 18]

