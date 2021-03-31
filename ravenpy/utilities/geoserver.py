"""
GeoServer interaction operations.

Working assumptions for this module:
* Point coordinates are passed as shapely.geometry.Point instances.
* BBox coordinates are passed as (lon1, lat1, lon2, lat2).
* Shapes (polygons) are passed as shapely.geometry.shape parsable objects.
* All functions that require a CRS have a CRS argument with a default set to WGS84.
* GEO_URL points to the GeoServer instance hosting all files.

"""
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple, Union
from urllib.parse import urljoin

from requests import Request

from . import gis_import_error_message

try:
    import fiona
    import pandas as pd
    from lxml import etree
    from owslib.fes import PropertyIsLike
    from owslib.wcs import WebCoverageService
    from owslib.wfs import WebFeatureService
    from shapely.geometry import Point, shape
except (ImportError, ModuleNotFoundError) as e:
    msg = gis_import_error_message.format(Path(__file__).stem)
    raise ImportError(msg) from e

# Do not remove the trailing / otherwise `urljoin` will remove the geoserver path.
GEO_URL = "http://pavics.ouranos.ca/geoserver/"

# We store the contour of different hydrobasins domains
hybas_dir = Path(__file__).parent.parent / "data" / "hydrobasins_domains"
hybas_pat = "hybas_lake_{}_lev01_v1c.zip"

# This could be inferred from existing files in hybas_dir
hybas_regions = ["na", "ar"]
hybas_domains = {dom: hybas_dir / hybas_pat.format(dom) for dom in hybas_regions}


def _get_location_wfs(
    coordinates: Tuple[
        Union[int, float, str],
        Union[str, float, int],
        Union[str, float, int],
        Union[str, float, int],
    ],
    layer: str = True,
    geoserver: str = GEO_URL,
) -> bytes:
    """Return leveled features from a hosted data set using bounding box coordinates and WFS 1.1.0 protocol.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    coordinates : Tuple[Union[str, float, int], Union[str, float, int], Union[str, float, int], Union[str, float, int]]
      Geographic coordinates of the bounding box (left, down, right, up).
    layer : str
      The WFS/WMS layer name requested.
    geoserver: str
      The address of the geoserver housing the layer to be queried. Default: http://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
      A GML-encoded vector feature.

    """
    wfs = WebFeatureService(url=urljoin(geoserver, "wfs"), version="1.1.0", timeout=30)
    resp = wfs.getfeature(
        typename=layer, bbox=coordinates, srsname="urn:x-ogc:def:crs:EPSG:4326"
    )

    data = resp.read()
    return data


def _get_feature_attributes_wfs(
    attribute: str = None,
    value: Union[str, float, int] = None,
    layer: str = None,
    geoserver: str = GEO_URL,
) -> str:
    """Return a URL that formats and returns remote GetFeatures request from a hosted WFS dataset.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    attribute : str
      Attribute/field to be queried.
    value: Union[str, float, int]
      Value for attribute queried.
    geoserver: str
      The address of the geoserver housing the layer to be queried. Default: http://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
      URL to the GeoJSON-encoded WFS response.

    """
    if attribute is None or value is None:
        raise NotImplementedError()

    try:
        attribute = str(attribute)
        value = str(value)

    except ValueError:
        raise Exception("Unable to cast attribute/filter to string")

    filter_request = PropertyIsLike(propertyname=attribute, literal=value, wildCard="*")
    filterxml = etree.tostring(filter_request.toXML()).decode("utf-8")
    params = dict(
        service="WFS",
        version="1.1.0",
        request="GetFeature",
        typename=layer,
        outputFormat="json",
        filter=filterxml,
    )

    q = Request("GET", url=urljoin(geoserver, "wfs"), params=params).prepare().url

    return q


def _determine_upstream_ids(
    fid: str,
    df: pd.DataFrame,
    basin_field: str = None,
    downstream_field: str = None,
    basin_family: Optional[str] = None,
) -> pd.Series:
    """Return a list of upstream features by evaluating the downstream networks.

    Parameters
    ----------
    fid : str
      feature ID of the downstream feature of interest.
    df : pd.DataFrame
      Dataframe comprising the watershed attributes.
    basin_field: str
      The field used to determine the id of the basin according to hydro project.
    downstream_field: str
      The field identifying the downstream sub-basin for the hydro project.
    basin_family: str, optional
      Regional watershed code (For HydroBASINS dataset).

    Returns
    -------
    pd.Series
      Basins ids including `fid` and its upstream contributors.
    """

    def upstream_ids(bdf, bid):
        return bdf[bdf[downstream_field] == bid][basin_field]

    # Note: Hydro Routing `SubId` is a float for some reason and Python float != GeoServer double. Cast them to int.
    if isinstance(fid, float):
        fid = int(fid)
        df[basin_field] = df[basin_field].astype(int)
        df[downstream_field] = df[downstream_field].astype(int)

    # Locate the downstream feature
    ds = df.set_index(basin_field).loc[fid]
    if basin_family is not None:
        # Do a first selection on the main basin ID of the downstream feature.
        sub = df[df[basin_family] == ds[basin_family]]
    else:
        sub = None

    # Find upstream basins
    up = [fid]
    for b in up:
        tmp = upstream_ids(sub if sub is not None else df, b)
        if len(tmp):
            up.extend(tmp)

    return (
        sub[sub[basin_field].isin(up)]
        if sub is not None
        else df[df[basin_field].isin(up)]
    )


def get_raster_wcs(
    coordinates: Union[Iterable, Sequence[Union[float, str]]],
    geographic: bool = True,
    layer: str = None,
    geoserver: str = GEO_URL,
) -> bytes:
    """Return a subset of a raster image from the local GeoServer via WCS 2.0.1 protocol.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    coordinates : Sequence[Union[int, float, str]]
      Geographic coordinates of the bounding box (left, down, right, up)
    geographic : bool
      If True, uses "Long" and "Lat" in WCS call. Otherwise uses "E" and "N".
    layer : str
      Layer name of raster exposed on GeoServer instance, e.g. 'public:CEC_NALCMS_LandUse_2010'
    geoserver: str
      The address of the geoserver housing the layer to be queried. Default: http://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    bytes
      A GeoTIFF array.

    """
    (left, down, right, up) = coordinates

    if geographic:
        x, y = "Long", "Lat"
    else:
        x, y = "E", "N"

    wcs = WebCoverageService(url=urljoin(geoserver, "ows"), version="2.0.1")

    try:
        resp = wcs.getCoverage(
            identifier=[layer],
            format="image/tiff",
            subsets=[(x, left, right), (y, down, up)],
            timeout=120,
        )

    except Exception as e:
        raise Exception(e)

    data = resp.read()

    try:
        etree.fromstring(data)
        # The response is an XML file describing the server error.
        raise ChildProcessError(data)

    except etree.XMLSyntaxError:
        # The response is the DEM array.
        return data


# ~~~~ HydroBASINS functions ~~~~ #


def hydrobasins_upstream_ids(
    fid: str,
    df: pd.DataFrame,
    basin_field="HYBAS_ID",
    downstream_field="NEXT_DOWN",
    basin_family="MAIN_BAS",
) -> pd.Series:
    """Return a list of hydrobasins features located upstream.

    Parameters
    ----------
    fid : str
      Basin feature ID code of the downstream feature.
    df : pd.DataFrame
      Watershed attributes.
    basin_field: str
      The field used to determine the id of the basin. Default: "HYBAS_ID".
    downstream_field: str
      The field identifying the downstream sub-basin. Default: "NEXT_DOWN".
    basin_family: str, optional
      Regional watershed code. Default: "MAIN_BAS".

    Returns
    -------
    pd.Series
      Basins ids including `fid` and its upstream contributors.
    """
    df_upstream = _determine_upstream_ids(
        fid=fid,
        df=df,
        basin_field=basin_field,
        downstream_field=downstream_field,
        basin_family=basin_family,
    )

    return df_upstream


def hydrobasins_aggregate(gdf: pd.DataFrame) -> pd.Series:
    """Aggregate multiple HydroBASINS watersheds into a single geometry.

    Parameters
    ----------
    gdf : pd.DataFrame
      Watershed attributes indexed by HYBAS_ID

    Returns
    -------
    pd.Series
    """

    # TODO: Review. Not sure it all makes sense. --> Looks fine to me? (TJS)
    def aggfunc(x):
        if x.name in ["COAST", "DIST_MAIN", "DIST_SINK"]:
            return x.min()
        elif x.name in ["SUB_AREA", "LAKE"]:
            return x.sum()
        else:
            return x[0]

    # Buffer function to fix invalid geometries
    gdf["geometry"] = gdf.buffer(0)

    return gdf.dissolve(by="MAIN_BAS", aggfunc=aggfunc)


def select_hybas_domain(
    bbox: Tuple[
        Union[int, float], Union[int, float], Union[int, float], Union[int, float]
    ]
) -> str:
    """
    Provided a given coordinate or boundary box, return the domain name of the geographic region
     the coordinate is located within.

    Parameters
    ----------
    bbox : tuple
      Geographic coordinates of the bounding box (left, down, right, up).

    Returns
    -------
    str
      The domain that the coordinate falls within. Possible results: "na", "ar".
    """

    for dom, fn in hybas_domains.items():
        with open(fn, "rb") as f:
            zf = fiona.io.ZipMemoryFile(f)
            coll = zf.open(fn.stem + ".shp")
            for _ in coll.filter(bbox=bbox):
                return dom

    raise LookupError(f"Could not find feature containing bbox: {bbox}.")


def get_hydrobasins_attributes_wfs(
    attribute: str = None,
    value: Union[str, float, int] = None,
    level: int = 12,
    lakes: bool = True,
    domain: str = None,
    geoserver: str = GEO_URL,
) -> str:
    """Return a URL that formats and returns a remote GetFeatures request from the USGS HydroBASINS dataset.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    attribute : str
      Attribute/field to be queried.
    value: Union[str, float, int]
      Value for attribute queried.
    level : int
      Level of granularity requested for the lakes vector (range(1,13)). Default: 12.
    lakes : bool
      Whether or not the vector should include the delimitation of lakes.
    domain : str
      The domain of the HydroBASINS data. Possible values:"na", "ar".
    geoserver: str
      The address of the geoserver housing the layer to be queried. Default: http://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
      URL to the GeoJSON-encoded WFS response.

    """
    layer = f"public:USGS_HydroBASINS_{'lake_' if lakes else ''}{domain}_lev{str(level).zfill(2)}"
    q = _get_feature_attributes_wfs(
        attribute=attribute, value=value, layer=layer, geoserver=geoserver
    )

    return q


def get_hydrobasins_location_wfs(
    coordinates: Tuple[
        Union[str, float, int],
        Union[str, float, int],
        Union[str, float, int],
        Union[str, float, int],
    ],
    level: int = 12,
    lakes: bool = True,
    domain: str = None,
    geoserver: str = GEO_URL,
) -> bytes:
    """Return features from the USGS HydroBASINS data set using bounding box coordinates.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    coordinates : Tuple[Union[str, float, int], Union[str, float, int], Union[str, float, int], Union[str, float, int]]
      Geographic coordinates of the bounding box (left, down, right, up).
    level : int
      Level of granularity requested for the lakes vector (1:12). Default: 12.
    lakes : bool
      Whether or not the vector should include the delimitation of lakes.
    domain : str
      The domain of the HydroBASINS data. Possible values:"na", "ar".
    geoserver: str
      The address of the geoserver housing the layer to be queried. Default: http://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
      A GML-encoded vector feature.

    """
    layer = f"public:USGS_HydroBASINS_{'lake_' if lakes else ''}{domain}_lev{str(level).zfill(2)}"
    data = _get_location_wfs(coordinates, layer=layer, geoserver=geoserver)

    return data


# ~~~~ Hydro Routing ~~~~ #


def hydro_routing_aggregate(gdf: pd.DataFrame) -> pd.Series:
    """Aggregate multiple hydro routing watersheds into a single geometry.

    Parameters
    ----------
    gdf : pd.DataFrame
      Watershed attributes indexed by HYBAS_ID

    Returns
    -------
    pd.Series
    """

    # TODO: @huard this needs to be discussed/reviewed. Dependent on ingesting the entire dataset. Very slow.
    def aggfunc(x):
        if x.name in ["area"]:
            return x.sum()
        elif x.name in ["MeanElev"]:
            return x.mean()
        else:
            return x[0]

    # Buffer function to fix invalid geometries
    gdf["geometry"] = gdf.buffer(0)

    return gdf.dissolve(aggfunc=aggfunc)


def hydro_routing_upstream_ids(
    fid: Union[str, float, int],
    df: pd.DataFrame,
    basin_field="SubId",
    downstream_field="DowSubId",
) -> pd.Series:
    """Return a list of hydro routing features located upstream.

    Parameters
    ----------
    fid : Union[str, float, int]
      Basin feature ID code of the downstream feature.
    df : pd.DataFrame
      Watershed attributes.
    basin_field: str
      The field used to determine the id of the basin. Default: "SubId".
    downstream_field: str
      The field identifying the downstream sub-basin. Default: "DowSubId".

    Returns
    -------
    pd.Series
      Basins ids including `fid` and its upstream contributors.
    """
    df_upstream = _determine_upstream_ids(
        fid=fid,
        df=df,
        basin_field=basin_field,
        downstream_field=downstream_field,
    )

    return df_upstream


def get_hydro_routing_attributes_wfs(
    attribute: str = None,
    value: Union[str, float, int] = None,
    level: int = 12,
    lakes: str = None,
    geoserver: str = GEO_URL,
) -> str:
    """Return a URL that formats and returns a remote GetFeatures request from hydro routing dataset.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    attribute : str
      Attribute/field to be queried.
    value: Union[str, float, int]
      Value for attribute queried.
    level : int
      Level of granularity requested for the lakes vector (range(7,13)). Default: 12.
    lakes : bool
      Query the version of dataset with lakes under 1km in width removed ("1km") or return all lakes ("all").
    geoserver: str
      The address of the geoserver housing the layer to be queried. Default: http://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
      URL to the GeoJSON-encoded WFS response.

    """
    layer = f"public:routing_{lakes}Lakes_{str(level).zfill(2)}"
    q = _get_feature_attributes_wfs(
        attribute=attribute, value=value, layer=layer, geoserver=geoserver
    )

    return q


def get_hydro_routing_location_wfs(
    coordinates: Tuple[
        Union[int, float, str],
        Union[str, float, int],
        Union[str, float, int],
        Union[str, float, int],
    ],
    lakes: str,
    level: int = 12,
    geoserver: str = GEO_URL,
) -> bytes:
    """Return features from the hydro routing data set using bounding box coordinates.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    coordinates : Tuple[Union[str, float, int], Union[str, float, int], Union[str, float, int], Union[str, float, int]]
      Geographic coordinates of the bounding box (left, down, right, up).
    lakes : {"1km", "all"}
      Query the version of dataset with lakes under 1km in width removed ("1km") or return all lakes ("all").
    level : int
      Level of granularity requested for the lakes vector (range(7,13)). Default: 12.
    geoserver: str
      The address of the geoserver housing the layer to be queried. Default: http://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
      A GML-encoded vector feature.

    """
    layer = f"public:routing_{lakes}Lakes_{str(level).zfill(2)}"
    data = _get_location_wfs(coordinates, layer=layer, geoserver=geoserver)

    return data
