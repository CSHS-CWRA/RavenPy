import tempfile

import numpy as np
import pytest

from ravenpy.utilities.testdata import get_local_testdata

pytestmark = pytest.mark.online

# FIXME: Remove XFAIL marks once OWSLib > 0.24.1 is released.


class TestHydroBASINS:
    geoserver = pytest.importorskip("ravenpy.utilities.geoserver")

    fiona = pytest.importorskip("fiona")
    gpd = pytest.importorskip("geopandas")
    sgeo = pytest.importorskip("shapely.geometry")

    def test_select_hybas_na_domain(self):
        bbox = (-68.0, 50.0) * 2
        dom = self.geoserver.select_hybas_domain(bbox)
        assert dom == "na"

    def test_select_hybas_ar_domain(self):
        bbox = (-114.65, 61.35) * 2
        dom = self.geoserver.select_hybas_domain(bbox)
        assert dom == "ar"

    def test_get_hydrobasins_location_wfs(self, tmp_path):
        lake_winnipeg = (
            -98.03575958286369,
            52.88238524279493,
        )
        resp = self.geoserver.get_hydrobasins_location_wfs(
            coordinates=lake_winnipeg, domain="na"
        )
        feat = self.gpd.GeoDataFrame.from_features(resp)
        geom = self.sgeo.shape(feat["geometry"][0])
        assert geom.bounds == (-99.2731, 50.3603, -96.2578, 53.8705)
        np.testing.assert_almost_equal(geom.area, 3.2530867)

    def test_get_hydrobasins_attributes_wfs(self, tmp_path):
        rio_grande = (-80.475, 8.4)
        resp = self.geoserver.get_hydrobasins_location_wfs(
            coordinates=rio_grande * 2, domain="na"
        )
        feat = self.gpd.GeoDataFrame.from_features(resp)
        main_bas = feat["MAIN_BAS"][0]

        region_url = self.geoserver.filter_hydrobasins_attributes_wfs(
            attribute="MAIN_BAS", value=main_bas, domain="na"
        )
        gdf = self.gpd.read_file(region_url)

        assert len(gdf) == 18
        assert gdf.crs.to_epsg() == 4326
        assert set(gdf["MAIN_BAS"].values) == {7120000210}

        basin = gdf.dissolve(by="MAIN_BAS")
        np.testing.assert_equal(
            basin.geometry.bounds.values,
            np.array([[-80.8542, 8.2459, -80.1375, 8.7004]]),
        )

    def test_hydrobasins_upstream_aggregate(self, tmp_path):
        puerto_cortes = (-83.525, 8.96, -83.520, 8.97)
        resp = self.geoserver.get_hydrobasins_location_wfs(
            coordinates=puerto_cortes, domain="na"
        )
        feat = self.gpd.GeoDataFrame.from_features(resp)

        gdf_upstream = self.geoserver.hydrobasins_upstream(feat.loc[0], domain="na")
        assert len(gdf_upstream) == 73
        aggregated = self.geoserver.hydrobasins_aggregate(gdf_upstream)

        assert len(aggregated) == 1
        assert float(aggregated.SUB_AREA.values) == 4977.8
        np.testing.assert_equal(
            aggregated.geometry.bounds.values,
            np.array([[-83.8167, 8.7625, -82.7125, 9.5875]]),
        )


class TestHydroRouting:
    geoserver = pytest.importorskip("ravenpy.utilities.geoserver")

    fiona = pytest.importorskip("fiona")
    gpd = pytest.importorskip("geopandas")
    sgeo = pytest.importorskip("shapely.geometry")

    def test_hydro_routing_locations(self, tmp_path):
        lake_winnipeg = (
            -98.03575958286369,
            52.88238524279493,
        )
        resp = self.geoserver.get_hydro_routing_location_wfs(
            coordinates=lake_winnipeg, lakes="all"
        )
        feat = self.gpd.GeoDataFrame.from_features(resp)
        geom = feat["geometry"][0]
        assert geom.bounds == (-99.3083, 50.1875, -95.9875, 54.0542)
        # Note: This value is not in sq. km.
        np.testing.assert_almost_equal(geom.area, 4.0978987)

    @pytest.mark.slow
    def test_get_hydro_routing_attributes_wfs(self):
        region_url = self.geoserver.filter_hydro_routing_attributes_wfs(
            attribute="IsLake", value="1.0", lakes="1km", level="07"
        )
        gdf = self.gpd.read_file(region_url)
        assert len(gdf) == 11415

    @pytest.mark.slow
    def test_hydro_routing_upstream(self, tmp_path):
        amadjuak = (-71.225, 65.05, -71.220, 65.10)
        resp = self.geoserver.get_hydro_routing_location_wfs(
            coordinates=amadjuak, lakes="1km", level=7
        )
        feat = self.gpd.GeoDataFrame.from_features(resp)
        subbasin_id = feat["SubId"][0]

        gdf_upstream = self.geoserver.hydro_routing_upstream(
            subbasin_id, lakes="1km", level=7
        )

        assert len(gdf_upstream) == 33  # TODO: Verify this with the model maintainers.


class TestWFS:
    geoserver = pytest.importorskip("ravenpy.utilities.geoserver")

    fiona = pytest.importorskip("fiona")
    gpd = pytest.importorskip("geopandas")
    sgeo = pytest.importorskip("shapely.geometry")

    def test_get_location_wfs_point(self, tmp_path):
        las_vegas = (-115.136389, 36.175)
        usa_admin_bounds = "public:usa_admin_boundaries"
        resp = self.geoserver._get_location_wfs(point=las_vegas, layer=usa_admin_bounds)
        feat = self.gpd.GeoDataFrame.from_features(resp)

        geom = feat["geometry"][0]
        assert geom.bounds == (-120.001, 35.0019, -114.0417, 41.9948)
        # Note: This value is not in sq. km.
        np.testing.assert_almost_equal(geom.area, 29.9690150)

    def test_get_location_wfs_bbox(self, tmp_path):
        new_vegas = (-115.697, 34.742, -114.279, 36.456)
        usa_admin_bounds = "public:usa_admin_boundaries"
        resp = self.geoserver._get_location_wfs(bbox=new_vegas, layer=usa_admin_bounds)
        feat = self.gpd.GeoDataFrame.from_features(resp)

        assert len(feat["geometry"]) == 3
        assert set(feat.STATE_NAME.unique()) == {"Nevada", "California", "Arizona"}

    def test_get_feature_attributes_wfs(self):
        state_name = "Nevada"
        usa_admin_bounds = "public:usa_admin_boundaries"

        vector_url = self.geoserver._filter_feature_attributes_wfs(
            attribute="STATE_NAME", value=state_name, layer=usa_admin_bounds
        )
        gdf = self.gpd.read_file(vector_url)
        assert len(gdf) == 1
        assert gdf.STATE_NAME.unique() == "Nevada"


class TestWCS:
    io = pytest.importorskip("ravenpy.utilities.io")
    geoserver = pytest.importorskip("ravenpy.utilities.geoserver")
    geo = pytest.importorskip("ravenpy.utilities.geo")
    rasterio = pytest.importorskip("rasterio")

    vector_file = get_local_testdata("polygons/Saskatoon.geojson")

    def test_get_raster_wcs(self, tmp_path):
        # TODO: This CRS needs to be redefined using modern pyproj-compatible strings.
        nalcms_crs = "+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs=True"

        with tempfile.NamedTemporaryFile(
            prefix="reprojected_", suffix=".json", dir=tmp_path
        ) as projected:
            self.geo.generic_vector_reproject(
                self.vector_file, projected.name, target_crs=nalcms_crs
            )
            bbox = self.io.get_bbox(projected.name)

        raster_url = "public:CEC_NALCMS_LandUse_2010"
        raster_bytes = self.geoserver.get_raster_wcs(
            bbox, geographic=False, layer=raster_url
        )

        with tempfile.NamedTemporaryFile(
            prefix="wcs_", suffix=".tiff", dir=tmp_path
        ) as raster_file:
            with open(raster_file.name, "wb") as rf:
                rf.write(raster_bytes)
                rf.close()
                with self.rasterio.open(rf.name) as src:
                    assert src.width == 650
                    assert src.height == 744
                    np.testing.assert_allclose(
                        src.lnglat(), (-106.641957, 52.156555), rtol=1e-5
                    )
                    data = src.read()
                    assert np.unique(data).tolist() == [1, 5, 8, 10, 14, 15, 16, 17, 18]
