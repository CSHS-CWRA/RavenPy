"""GeoServer interaction operations.

Working assumptions for this module:
* Point coordinates are passed as shapely.geometry.Point instances.
* BBox coordinates are passed as (lon1, lat1, lon2, lat2).
* Shapes (polygons) are passed as shapely.geometry.shape parsable objects.
* All functions that require a CRS have a CRS argument with a default set to WGS84.
* GEOSERVER_URL points to the GeoServer instance hosting all files.
* For legacy reasons, we also accept the `GEO_URL` environment variable.

TODO: Refactor to remove functions that are just 2-lines of code.
For example, many function's logic essentially consists in creating the layer name.
We could have a function that returns the layer name, and then other functions expect the layer name.
"""

import inspect
import json
import os
import warnings
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Optional, Union
from urllib.parse import urlencode, urljoin

import urllib3

from . import gis_import_error_message

try:
    import fiona
    import geopandas as gpd
    import pandas as pd
    from lxml import etree
    from owslib.fes import PropertyIsLike
    from owslib.wcs import WebCoverageService
    from owslib.wfs import WebFeatureService
    from shapely.geometry import Point, shape
except (ImportError, ModuleNotFoundError) as e:
    msg = gis_import_error_message.format(Path(__file__).stem)
    raise ImportError(msg) from e

try:
    from owslib.fes2 import Intersects
    from owslib.gml import Point as wfs_Point
except (ImportError, ModuleNotFoundError):
    warnings.warn("WFS point spatial filtering requires OWSLib>0.24.1.")
    Intersects = None
    wfs_Point = None  # noqa: N816

from .geo import determine_upstream_ids

# Can be set at runtime with `$ env RAVENPY_GEOSERVER_URL=https://xx.yy.zz/geoserver/ ...`.
# For legacy reasons, we also accept the `GEO_URL` environment variable.
HOST_URL = os.getenv(
    "RAVENPY_HOST_URL", os.getenv("HOST_URL", "https://pavics.ouranos.ca/")
)
GEOSERVER_URL = os.getenv(
    "RAVENPY_GEOSERVER_URL",
    os.getenv("GEO_URL", f"{HOST_URL}/geoserver/"),
)
if not GEOSERVER_URL.endswith("/"):
    GEOSERVER_URL = f"{GEOSERVER_URL}/"

# We store the contour of different HydroBASINS domains
hybas_dir = Path(__file__).parent.parent / "data" / "hydrobasins_domains"
hybas_pat = "hybas_lake_{domain}_lev01_v1c.zip"

# This could be inferred from existing files in hybas_dir
hybas_regions = ["na", "ar"]
hybas_domains = {dom: hybas_dir / hybas_pat.format(domain=dom) for dom in hybas_regions}


def _fix_server_url(server_url: str) -> str:
    if not server_url.endswith("/"):
        warnings.warn("The url should end with a slash. Appending it to the url.")
        return f"{server_url}/"
    return server_url


def _get_location_wfs(
    bbox: Optional[
        tuple[
            Union[str, float, int],
            Union[str, float, int],
            Union[str, float, int],
            Union[str, float, int],
        ]
    ] = None,
    point: Optional[
        tuple[
            Union[str, float, int],
            Union[str, float, int],
        ]
    ] = None,
    *,
    layer: str,
    geoserver: str = GEOSERVER_URL,
) -> dict:
    """Return leveled features from a hosted data set using bounding box coordinates and WFS 1.1.0 protocol.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    bbox : Optional[Tuple[Union[str, float, int], Union[str, float, int], Union[str, float, int], Union[str, float, int]]]
        Geographic coordinates of the bounding box (left, down, right, up).
    point : Optional[Tuple[Union[str, float, int], Union[str, float, int]]]
        Geographic coordinates of an intersecting point (lon, lat).
    layer : str
        The WFS/WMS layer name requested.
    geoserver : str
        The address of the geoserver housing the layer to be queried. Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    dict
        A GeoJSON-derived dictionary of vector features (FeatureCollection).
    """
    geoserver = _fix_server_url(geoserver)

    wfs = WebFeatureService(url=urljoin(geoserver, "wfs"), version="2.0.0", timeout=30)

    if bbox and point:
        raise NotImplementedError("Provide either 'bbox' or 'point'.")
    if bbox:
        kwargs = dict(bbox=bbox)
    elif point:
        # FIXME: Remove this once OWSlib > 0.24.1 is released.
        if not Intersects and not wfs_Point:
            raise NotImplementedError(
                f"{inspect.stack()[1][3]} with point filtering requires OWSLib>0.24.1.",
            )

        p = wfs_Point(
            id="feature",
            srsName="http://www.opengis.net/gml/srs/epsg.xml#4326",
            pos=point,
        )
        f = Intersects(propertyname="the_geom", geometry=p)
        intersects = f.toXML()
        kwargs = dict(filter=intersects)
    else:
        raise ValueError()

    resp = wfs.getfeature(
        typename=layer, outputFormat="application/json", method="POST", **kwargs
    )

    data = json.loads(resp.read())
    return data


def _get_feature_attributes_wfs(
    attribute: Sequence[str],
    layer: str,
    geoserver: str = GEOSERVER_URL,
) -> str:
    """Return WFS GetFeature URL request for attribute values.

    Making this request will return a JSON response.

    Parameters
    ----------
    attribute : list
        Attribute/field names.
    layer : str
        Name of geographic layer queried.
    geoserver : str
        The address of the geoserver housing the layer to be queried.
        Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
        WFS request URL.

    Notes
    -----
    Non-existent attributes will raise a cryptic DriverError from fiona.
    """
    host = _fix_server_url(HOST_URL)
    geoserver = _fix_server_url(geoserver)

    params = dict(
        service="WFS",
        version="2.0.0",
        request="GetFeature",
        typename=layer,
        outputFormat="application/json",
        propertyName=",".join(attribute),
    )
    url = urljoin(geoserver, "wfs") + "?" + urlencode(params)
    http = urllib3.PoolManager()
    response = http.request("GET", url)

    return urljoin(host, response.url)


def _filter_feature_attributes_wfs(
    attribute: str,
    value: Union[str, float, int],
    layer: str,
    geoserver: str = GEOSERVER_URL,
) -> str:
    """Return WFS GetFeature URL request filtering geographic features based on a property's value.

    Parameters
    ----------
    attribute : str
        Attribute/field name.
    value: Union[str, float, int]
        Value for attribute queried.
    layer : str
        Name of geographic layer queried.
    geoserver : str
        The address of the geoserver housing the layer to be queried.
        Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
      WFS request URL.
    """
    host = _fix_server_url(HOST_URL)
    geoserver = _fix_server_url(geoserver)

    try:
        attribute = str(attribute)
        value = str(value)

    except ValueError:
        raise Exception("Unable to cast attribute/filter to string.")

    filter_request = PropertyIsLike(propertyname=attribute, literal=value, wildCard="*")
    filterxml = etree.tostring(filter_request.toXML()).decode("utf-8")
    params = dict(
        service="WFS",
        version="1.1.0",
        request="GetFeature",
        typename=layer,
        outputFormat="application/json",
        filter=filterxml,
    )

    url = urljoin(geoserver, "wfs") + "?" + urlencode(params)
    http = urllib3.PoolManager()
    response = http.request("GET", url)

    return urljoin(host, response.url)


def get_raster_wcs(
    coordinates: Union[Iterable, Sequence[Union[float, str]]],
    geographic: bool = True,
    *,
    layer: str,
    geoserver: str = GEOSERVER_URL,
) -> bytes:
    """Return a subset of a raster image from the local GeoServer via WCS 2.0.1 protocol.

    For geographic raster grids, subsetting is based on WGS84 (Long, Lat) boundaries.
    If not geographic, subsetting based on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    coordinates : Sequence of int or float or str
        Geographic coordinates of the bounding box (left, down, right, up).
    geographic : bool
        If True, uses "Long" and "Lat" in WCS call. Otherwise, uses "E" and "N".
    layer : str
        Layer name of raster exposed on GeoServer instance, e.g. 'public:CEC_NALCMS_LandUse_2010'.
    geoserver : str
        The address of the geoserver housing the layer to be queried.
        Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    bytes
        A GeoTIFF array.
    """
    geoserver = _fix_server_url(geoserver)

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

    except Exception:
        raise

    data = resp.read()

    try:
        etree.fromstring(data)  # noqa: S320
        # The response is an XML file describing the server error.
        raise ChildProcessError(data)

    except etree.XMLSyntaxError:
        # The response is the DEM array.
        return data


# ~~~~ HydroBASINS functions ~~~~ #


def hydrobasins_upstream(feature: dict, domain: str) -> pd.DataFrame:
    """Return a list of HydroBASINS features located upstream.

    Parameters
    ----------
    feature : dict
        Basin feature attributes, including the fields ["HYBAS_ID", "NEXT_DOWN", "MAIN_BAS"].
    domain : {"na", "ar"}
        Domain of the feature, North America or Arctic.

    Returns
    -------
    pd.Series
        Basins ids including `fid` and its upstream contributors.
    """
    basin_field = "HYBAS_ID"
    downstream_field = "NEXT_DOWN"
    basin_family = "MAIN_BAS"

    # This does not work with `wfs.getfeature`. No filtering occurs when asking for specific attributes.
    # wfs = WebFeatureService(url=urljoin(geoserver, "wfs"), version="2.0.0", timeout=30)
    # layer = f"public:USGS_HydroBASINS_{'lake_' if lakes else ''}{domain}_lev{str(level).zfill(2)}"
    # filter = PropertyIsEqualTo(propertyname=basin_family, literal=feature[basin_family])

    # Fetch all features in the same basin
    request_url = filter_hydrobasins_attributes_wfs(
        attribute=basin_family, value=feature[basin_family], domain=domain
    )
    df = gpd.read_file(filename=request_url, engine="pyogrio")

    # Filter upstream watersheds
    return determine_upstream_ids(
        fid=feature[basin_field],
        df=df,
        basin_field=basin_field,
        downstream_field=downstream_field,
    )


def hydrobasins_aggregate(gdf: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate multiple HydroBASINS watersheds into a single geometry.

    Parameters
    ----------
    gdf : pd.DataFrame
        Watershed attributes indexed by HYBAS_ID.

    Returns
    -------
    pd.DataFrame
    """
    i0 = gdf.index[0]

    # TODO: Review. Not sure it all makes sense. --> Looks fine to me? (TJS)
    def aggfunc(x):
        if x.name in ["COAST", "DIST_MAIN", "DIST_SINK"]:
            return x.min()
        elif x.name in ["SUB_AREA", "LAKE"]:
            return x.sum()
        else:
            return x.loc[i0]

    # Buffer function to fix invalid geometries
    gdf["geometry"] = gdf.buffer(0)

    return gdf.dissolve(by="MAIN_BAS", aggfunc=aggfunc)


def select_hybas_domain(
    bbox: Optional[
        tuple[
            Union[int, float], Union[int, float], Union[int, float], Union[int, float]
        ]
    ] = None,
    point: Optional[tuple[Union[int, float], Union[int, float]]] = None,
) -> str:
    """
    Provided a given coordinate or boundary box, return the domain name of the geographic region the coordinate is located within.

    Parameters
    ----------
    bbox : Optional[Tuple[Union[float, int], Union[float, int], Union[float, int], Union[float, int]]]
        Geographic coordinates of the bounding box (left, down, right, up).
    point : Optional[Tuple[Union[float, int], Union[float, int]]]
        Geographic coordinates of an intersecting point (lon, lat).

    Returns
    -------
    str
        The domain that the coordinate falls within. Possible results: "na", "ar".
    """
    if bbox and point:
        raise NotImplementedError("Provide either 'bbox' or 'point'.")
    if point:
        bbox = point * 2

    for dom, fn in hybas_domains.items():
        with Path(fn).open("rb") as f:
            zf = fiona.io.ZipMemoryFile(f)
            coll = zf.open(fn.stem + ".shp")
            for feat in coll.filter(bbox=bbox):
                if isinstance(feat, fiona.Feature):
                    return dom

    raise LookupError(f"Could not find feature containing bbox: {bbox}.")


def filter_hydrobasins_attributes_wfs(
    attribute: str,
    value: Union[str, float, int],
    domain: str,
    geoserver: str = GEOSERVER_URL,
) -> str:
    """
    Return a URL that formats and returns a remote GetFeatures request from the USGS HydroBASINS dataset.

    For geographic raster grids, subsetting is based on WGS84 (Long, Lat) boundaries.
    If not geographic, subsetting based on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    attribute : str
        Attribute/field to be queried.
    value : str or float or int
        Value for attribute queried.
    domain : {"na", "ar"}
        The domain of the HydroBASINS data.
    geoserver : str
        The address of the geoserver housing the layer to be queried. Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
        URL to the GeoJSON-encoded WFS response.
    """
    geoserver = _fix_server_url(geoserver)

    lakes = True
    level = 12

    layer = f"public:USGS_HydroBASINS_{'lake_' if lakes else ''}{domain}_lev{str(level).zfill(2)}"
    q = _filter_feature_attributes_wfs(
        attribute=attribute, value=value, layer=layer, geoserver=geoserver
    )

    return q


def get_hydrobasins_location_wfs(
    coordinates: tuple[
        Union[str, float, int],
        Union[str, float, int],
    ],
    domain: str,
    geoserver: str = GEOSERVER_URL,
) -> dict[str, Union[str, int, float]]:
    """
    Return features from the USGS HydroBASINS data set using bounding box coordinates.

    For geographic raster grids, subsetting is based on WGS84 (Long, Lat) boundaries.
    If not geographic, subsetting based on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    coordinates : Tuple[str or float or int, str or float or int]
        Geographic coordinates of the bounding box (left, down, right, up).
    domain : {"na", "ar"}
        The domain of the HydroBASINS data.
    geoserver : str
        The address of the geoserver housing the layer to be queried. Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    dict
        A GeoJSON-encoded vector feature.
    """
    geoserver = _fix_server_url(geoserver)

    lakes = True
    level = 12
    layer = f"public:USGS_HydroBASINS_{'lake_' if lakes else ''}{domain}_lev{str(level).zfill(2)}"
    if not wfs_Point and not Intersects:
        data = _get_location_wfs(bbox=coordinates * 2, layer=layer, geoserver=geoserver)
    else:
        data = _get_location_wfs(point=coordinates, layer=layer, geoserver=geoserver)

    return data


# ~~~~ Hydro Routing ~~~~ #


def hydro_routing_upstream(
    fid: Union[str, float, int],
    level: int = 12,
    lakes: str = "1km",
    geoserver: str = GEOSERVER_URL,
) -> pd.Series:
    """
    Return a list of hydro routing features located upstream.

    Parameters
    ----------
    fid : str or float or int
        Basin feature ID code of the downstream feature.
    level : int
        Level of granularity requested for the lakes vector (range(7,13)). Default: 12.
    lakes : {"1km", "all"}
        Query the version of dataset with lakes under 1km in width removed ("1km") or return all lakes ("all").
    geoserver : str
        The address of the geoserver housing the layer to be queried.
        Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    pd.Series
        Basins ids including `fid` and its upstream contributors.
    """
    geoserver = _fix_server_url(geoserver)

    wfs = WebFeatureService(url=urljoin(geoserver, "wfs"), version="2.0.0", timeout=30)
    layer = f"public:routing_{lakes}Lakes_{str(level).zfill(2)}"

    # Get attributes necessary to identify upstream watersheds
    resp = wfs.getfeature(
        typename=layer,
        propertyname=["SubId", "DowSubId"],
        outputFormat="application/json",
    )
    df = gpd.read_file(resp)

    # Identify upstream features
    df_upstream = determine_upstream_ids(
        fid=fid,
        df=df,
        basin_field="SubId",
        downstream_field="DowSubId",
    )

    # Fetch upstream features
    resp = wfs.getfeature(
        typename=layer,
        featureid=df_upstream["id"].tolist(),
        outputFormat="application/json",
    )

    return gpd.read_file(resp.read().decode())


def get_hydro_routing_attributes_wfs(
    attribute: Sequence[str],
    level: int = 12,
    lakes: str = "1km",
    geoserver: str = GEOSERVER_URL,
) -> str:
    """
    Return a URL that formats and returns a remote GetFeatures request from hydro routing dataset.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries.
    If not geographic, subsetting based on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    attribute : list
        Attributes/fields to be queried.
    level : int
        Level of granularity requested for the lakes vector (range(7,13)). Default: 12.
    lakes : {"1km", "all"}
        Query the version of dataset with lakes under 1km in width removed ("1km") or return all lakes ("all").
    geoserver : str
        The address of the geoserver housing the layer to be queried.
        Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
        URL to the GeoJSON-encoded WFS response.
    """
    geoserver = _fix_server_url(geoserver)

    layer = f"public:routing_{lakes}Lakes_{str(level).zfill(2)}"
    return _get_feature_attributes_wfs(
        attribute=attribute, layer=layer, geoserver=geoserver
    )


def filter_hydro_routing_attributes_wfs(
    attribute: str,
    value: Union[str, float, int],
    level: int = 12,
    lakes: str = "1km",
    geoserver: str = GEOSERVER_URL,
) -> str:
    """
    Return a URL that formats and returns a remote GetFeatures request from hydro routing dataset.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    attribute : list
        Attributes/fields to be queried.
    value : str or int or float
        The requested value for the attribute.
    level : int
        Level of granularity requested for the lakes vector (range(7,13)). Default: 12.
    lakes : {"1km", "all"}
        Query the version of dataset with lakes under 1km in width removed ("1km") or return all lakes ("all").
    geoserver : str
        The address of the geoserver housing the layer to be queried.
        Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    str
        URL to the GeoJSON-encoded WFS response.
    """
    geoserver = _fix_server_url(geoserver)

    layer = f"public:routing_{lakes}Lakes_{str(level).zfill(2)}"
    return _filter_feature_attributes_wfs(
        attribute=attribute, value=value, layer=layer, geoserver=geoserver
    )


def get_hydro_routing_location_wfs(
    coordinates: tuple[
        Union[int, float, str],
        Union[str, float, int],
    ],
    lakes: str,
    level: int = 12,
    geoserver: str = GEOSERVER_URL,
) -> dict:
    """
    Return features from the hydro routing data set using bounding box coordinates.

    For geographic rasters, subsetting is based on WGS84 (Long, Lat) boundaries. If not geographic, subsetting based
    on projected coordinate system (Easting, Northing) boundaries.

    Parameters
    ----------
    coordinates : Tuple[str or float or int, str or float or int]
        Geographic coordinates of the bounding box (left, down, right, up).
    lakes : {"1km", "all"}
        Query the version of dataset with lakes under 1km in width removed ("1km") or return all lakes ("all").
    level : int
        Level of granularity requested for the lakes vector (range(7,13)). Default: 12.
    geoserver : str
        The address of the geoserver housing the layer to be queried.
        Default: https://pavics.ouranos.ca/geoserver/.

    Returns
    -------
    dict
        A GeoJSON-derived dictionary of vector features (FeatureCollection).
    """
    geoserver = _fix_server_url(geoserver)

    layer = f"public:routing_{lakes}Lakes_{str(level).zfill(2)}"

    if not wfs_Point and not Intersects:
        data = _get_location_wfs(bbox=coordinates * 2, layer=layer, geoserver=geoserver)
    else:
        data = _get_location_wfs(point=coordinates, layer=layer, geoserver=geoserver)

    return data
