import collections
import json
import logging
import math
import re
import tarfile
import tempfile
import warnings
import zipfile
from pathlib import Path
from re import search
from typing import Any, Iterable, List, Optional, Sequence, Tuple, Union

try:
    import fiona
    import fiona.crs
    import pyproj
    import rasterio
    import rasterio.mask
    import rasterio.vrt
    import rasterio.warp
    from affine import Affine
    from osgeo.gdal import DEMProcessing, Dataset
    from osgeo import gdal_array
    from pyproj.crs import CRS, CRSError
    from shapely.geometry import (
        GeometryCollection,
        MultiPolygon,
        Point,
        Polygon,
        mapping,
        shape,
    )
    from shapely.ops import transform
except (ImportError, ModuleNotFoundError) as e:
    msg = (
        f"`{Path(__file__).stem}` requires installation of the RavenPy GIS libraries. These can be installed using the"
        " `pip install ravenpy[gis]` recipe or via Anaconda (`conda env -n ravenpy-env -f environment.yml`)"
        " from the RavenPy repository source files."
    )
    raise ImportError(msg) from e

import numpy as np

LOGGER = logging.getLogger("RavenPy")

# See: https://kokoalberti.com/articles/geotiff-compression-optimization-guide/
GDAL_TIFF_COMPRESSION_OPTION = "compress=lzw"  # or 'compress=deflate' or 'compress=zstd' or 'compress=lerc' or others
RASTERIO_TIFF_COMPRESSION = "lzw"

WGS84 = 4326


def address_append(address: Union[str, Path]) -> str:
    """
    Formats a URL/URI to be more easily read with libraries such as "rasterstats"

    Parameters
    ----------
    address: Union[str, Path]
      URL/URI to a potential zip or tar file

    Returns
    -------
    str
      URL/URI prefixed for archive type
    """
    zipped = search(r"(\.zip)", str(address))
    tarred = search(r"(\.tar)", str(address))

    try:
        if zipped:
            return f"zip://{address}"
        elif tarred:
            return f"tar://{address}"
        else:
            LOGGER.info("No changes made to address.")
            return str(address)
    except Exception:
        LOGGER.error("Failed to prefix or parse URL %s." % address)


def generic_extract_archive(
    resources: Union[str, Path, List[Union[bytes, str, Path]]],
    output_dir: Optional[Union[str, Path]] = None,
) -> List[str]:
    """Extract archives (tar/zip) to a working directory.

    Parameters
    ----------
    resources: Union[str, Path, List[Union[bytes, str, Path]]]
      list of archive files (if netCDF files are in list, they are passed and returned as well in the return).
    output_dir: Optional[Union[str, Path]]
      string or Path to a working location (default: temporary folder).

    Returns
    -------
    list
      List of original or of extracted files
    """

    archive_types = [".tar", ".zip", ".7z"]
    output_dir = output_dir or tempfile.gettempdir()

    if not isinstance(resources, list):
        resources = [resources]

    files = list()

    for arch in resources:
        if any(ext in str(arch).lower() for ext in archive_types):
            try:
                LOGGER.debug("archive=%s", arch)
                file = Path(arch).name

                if file.endswith(".nc"):
                    files.append(Path(output_dir.join(arch)))
                elif file.endswith(".tar"):
                    with tarfile.open(arch, mode="r") as tar:
                        tar.extractall(path=output_dir)
                        files.extend(
                            [str(Path(output_dir).joinpath(f)) for f in tar.getnames()]
                        )
                elif file.endswith(".zip"):
                    with zipfile.ZipFile(arch, mode="r") as zf:
                        zf.extractall(path=output_dir)
                        files.extend(
                            [str(Path(output_dir).joinpath(f)) for f in zf.namelist()]
                        )
                elif file.endswith(".7z"):
                    msg = "7z file extraction is not supported at this time."
                    LOGGER.warning(msg)
                    warnings.warn(msg, UserWarning)
                else:
                    LOGGER.debug('File extension "%s" unknown' % file)
            except Exception as e:
                LOGGER.error("Failed to extract sub archive {%s}: {%s}" % (arch, e))
        else:
            LOGGER.warning("No archives found. Continuing...")
            return resources

    return files


def archive_sniffer(
    archives: Union[str, Path, List[Union[str, Path]]],
    working_dir: Optional[Union[str, Path]] = None,
    extensions: Optional[Iterable[str]] = None,
) -> List[Union[str, Path]]:
    """Return a list of locally unarchived files that match the desired extensions.

    Parameters
    ----------
    archives : Union[str, Path, List[Union[str, Path]]]
      archive location or list of archive locations
    working_dir : Union[str, path]
      string or Path to a working location
    extensions : List[str]
      list of accepted extensions

    Returns
    -------
    List[Union[str, Path]]
      List of files with matching accepted extensions
    """
    potential_files = list()

    if not extensions:
        extensions = [".gml", ".shp", ".geojson", ".gpkg", ".json"]

    decompressed_files = generic_extract_archive(archives, output_dir=working_dir)
    for file in decompressed_files:
        if any(ext in Path(file).suffix for ext in extensions):
            potential_files.append(file)
    return potential_files


def crs_sniffer(
    *args: Union[str, Path, Sequence[Union[str, Path]]]
) -> Union[List[Union[str, int]], str, int]:
    """Return the list of CRS found in files.

    Parameters
    ----------
    args : Union[str, Path, Sequence[Union[str, Path]]]
      Path(s) to the file(s) to examine.

    Returns
    -------
    Union[List[str], str]
      Returns either a list of CRSes or a single CRS definition, depending on the number of instances found.
    """
    crs_list = list()
    vectors = (".gml", ".shp", ".geojson", ".gpkg", ".json")
    rasters = (".tif", ".tiff")
    all_files = vectors + rasters

    for file in args:
        found_crs = False
        suffix = Path(file).suffix.lower()
        try:
            if suffix == ".zip":
                file = archive_sniffer(file, extensions=all_files)[0]
                suffix = Path(file).suffix.lower()

            if suffix in vectors:
                if suffix == ".gpkg":
                    if len(fiona.listlayers(file)) > 1:
                        raise NotImplementedError
                with fiona.open(file, "r") as src:
                    found_crs = CRS.from_wkt(src.crs_wkt).to_epsg()
            elif suffix in rasters:
                with rasterio.open(file, "r") as src:
                    found_crs = CRS.from_user_input(src.crs).to_epsg()
            else:
                raise FileNotFoundError("Invalid filename suffix")
        except FileNotFoundError as e:
            msg = f"{e}: Unable to open file {args}"
            LOGGER.warning(msg)
            raise Exception(msg)
        except NotImplementedError as e:
            msg = f"{e}: Multilayer GeoPackages are currently unsupported"
            LOGGER.error(msg)
            raise Exception(msg)
        except RuntimeError:
            pass

        crs_list.append(found_crs)

    if crs_list is None:
        msg = f"No CRS definitions found in {args}."
        raise FileNotFoundError(msg)

    if len(crs_list) == 1:
        if not crs_list[0]:
            msg = f"No CRS definitions found in {args}. Assuming {WGS84}."
            LOGGER.warning(msg)
            warnings.warn(msg, UserWarning)
            return WGS84
        return crs_list[0]
    return crs_list


def raster_datatype_sniffer(file: Union[str, Path]) -> str:
    """Return the type of the raster stored in the file.

    Parameters
    ----------
    file : Union[str, Path]
      Path to file.

    Returns
    -------
    str
      rasterio datatype of array values
    """
    try:
        with rasterio.open(file, "r") as src:
            dtype = src.dtypes[0]
        return dtype
    except rasterio.errors.RasterioError:
        msg = "Unable to read data type from {}.".format(file)
        LOGGER.exception(msg)
        raise ValueError(msg)


def parse_lonlat(lonlat: Union[str, Tuple[str, str]]) -> Tuple[float, float]:
    """Return longitude and latitude from a string.

    Parameters
    ----------
    lonlat : Union[str, Tuple[str, str]]
      A tuple or a str of lon and lat coordinates.

    Returns
    -------
    Tuple[float, float]
    """
    try:
        if isinstance(lonlat, str):
            lon, lat = tuple(map(float, re.findall(r"[-+]?[0-9]*\.?[0-9]+", lonlat)))
        elif isinstance(lonlat, tuple):
            lon, lat = map(float, lonlat)
        else:
            raise ValueError
        return lon, lat
    except Exception as e:
        msg = "Failed to parse longitude, latitude coordinates {}".format(lonlat)
        raise Exception(msg) from e


def single_file_check(file_list: List[Union[str, Path]]) -> Any:
    """Return the first element of a file list. Raise an error if the list is empty or contains more than one element.

    Parameters
    ----------
    file_list : List[Union[str, Path]]
    """
    if isinstance(file_list, (str, Path)):
        return file_list

    try:
        if len(file_list) > 1:
            msg = "Multi-file handling for file is not supported. Exiting."
            raise NotImplementedError(msg)
        elif len(file_list) == 0:
            msg = "No files found. Exiting."
            raise FileNotFoundError(msg)
        return file_list[0]
    except (FileNotFoundError, NotImplementedError) as e:
        LOGGER.error(e)
        raise


def boundary_check(
    *args: Sequence[Union[str, Path]],
    max_y: Union[int, float] = 60,
    min_y: Union[int, float] = -60,
) -> None:
    """Verify that boundaries do not exceed specific latitudes for geographic coordinate data. Raise a warning if so.

    Parameters
    ----------
    *args : Sequence[Union[str, Path]]
      listing of strings or paths to files
    max_y : Union[int, float]
      Maximum value allowed for latitude. Default: 60.
    min_y : Union[int, float]
      Minimum value allowed for latitude. Default: -60.
    """
    vectors = (".gml", ".shp", ".geojson", ".gpkg", ".json")
    rasters = (".tif", ".tiff")

    if len(args) == 1 and not isinstance(args[0], str):
        args = args[0]

    for file in args:
        try:
            if str(file).lower().endswith(vectors):
                src = fiona.open(file, "r")
            elif str(file).lower().endswith(rasters):
                src = rasterio.open(file, "r")
            else:
                raise FileNotFoundError()

            try:
                geographic = CRS(src.crs).is_geographic
            except CRSError:
                geographic = True
            src_min_y, src_max_y = src.bounds[1], src.bounds[3]
            if geographic and (src_max_y > max_y or src_min_y < min_y):
                msg = (
                    f"Vector {file} contains geometries in high latitudes."
                    " Verify choice of projected CRS is appropriate for analysis."
                )
                LOGGER.warning(msg)
                warnings.warn(msg, UserWarning)
            if not geographic:
                msg = f"Vector {file} is not in a geographic coordinate system."
                LOGGER.warning(msg)
                warnings.warn(msg, UserWarning)
            src.close()

        except FileNotFoundError:
            msg = f"Unable to read boundaries from {file}"
            LOGGER.error(msg)
            raise
    return


def multipolygon_check(geom: GeometryCollection) -> None:
    """Perform a check to verify a geometry is a MultiPolygon

    Parameters
    ----------
    geom : GeometryCollection

    Returns
    -------
    None
    """
    if not isinstance(type(geom), GeometryCollection):
        try:
            geom = shape(geom)
        except AttributeError:
            LOGGER.error("Unable to load argument as shapely.geometry.shape().")
            raise

    if isinstance(type(geom), MultiPolygon):
        LOGGER.warning("Shape is a Multipolygon.")
    return


def geom_transform(
    geom: Union[GeometryCollection, shape],
    source_crs: Union[str, int, CRS] = WGS84,
    target_crs: Union[str, int, CRS] = None,
) -> GeometryCollection:
    """Change the projection of a geometry.

    Assuming a geometry's coordinates are in a `source_crs`, compute the new coordinates under the `target_crs`.

    Parameters
    ----------
    geom : Union[GeometryCollection, shape]
      Source geometry.
    source_crs : Union[str, int, CRS]
      Projection identifier (proj4) for the source geometry, e.g. '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, int, CRS]
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    GeometryCollection
      Reprojected geometry.
    """
    try:
        from pyproj import Transformer  # noqa
        from functools import partial

        source = (
            CRS.from_epsg(source_crs)
            if isinstance(source_crs, int or str)
            else source_crs
        )
        target = (
            CRS.from_epsg(target_crs)
            if isinstance(target_crs, int or str)
            else target_crs
        )

        transform_func = Transformer.from_crs(source, target, always_xy=True)
        reprojected = transform(transform_func.transform, geom)

        return reprojected
    except Exception as err:
        msg = f"{err}: Failed to reproject geometry"
        LOGGER.error(msg)
        raise Exception(msg)


def geom_prop(geom: Union[Polygon, MultiPolygon, GeometryCollection]) -> dict:
    """Return a dictionary of geometry properties.

    Parameters
    ----------
    geom : Union[Polygon, MultiPolygon, GeometryCollection]
      Geometry to analyze.

    Returns
    -------
    dict
      Dictionary storing polygon area, centroid location, perimeter and gravelius shape index.

    Notes
    -----
    Some of the properties should be computed using an equal-area projection.
    """

    geom = shape(geom)
    lon, lat = geom.centroid.x, geom.centroid.y
    if (lon > 180) or (lon < -180) or (lat > 90) or (lat < -90):
        LOGGER.warning("Shape centroid is not in decimal degrees.")
    area = geom.area
    length = geom.length
    gravelius = length / 2 / math.sqrt(math.pi * area)
    parameters = {
        "area": area,
        "centroid": (lon, lat),
        "perimeter": length,
        "gravelius": gravelius,
    }
    return parameters


def dem_prop(
    dem: Union[str, Path],
    geom: Union[Polygon, MultiPolygon, List[Union[Polygon, MultiPolygon]]] = None,
    directory: Union[str, Path] = None,
) -> dict:
    """Return raster properties for each geometry.

    This

    Parameters
    ----------
    dem : Union[str, Path]
      DEM raster in reprojected coordinates.
    geom : Union[Polygon, MultiPolygon, List[Union[Polygon, MultiPolygon]]]
      Geometry over which aggregate properties will be computed. If None compute properties over entire raster.
    directory : Union[str, Path]
      Folder to save the GDAL terrain analysis outputs.

    Returns
    -------
    dict
      Dictionary storing mean elevation [m], slope [deg] and aspect [deg].
    """

    fns = dict()
    fns["dem"] = (
        tempfile.NamedTemporaryFile(
            prefix="dem", suffix=".tiff", dir=directory, delete=False
        ).name
        if geom is not None
        else dem
    )
    for key in ["slope", "aspect"]:
        fns[key] = tempfile.NamedTemporaryFile(
            prefix=key, suffix=".tiff", dir=directory, delete=False
        ).name

    # Clip to relevant area or read original raster
    if geom is None:
        with rasterio.open(dem) as f:
            elevation = f.read(1, masked=True)
    else:
        generic_raster_clip(raster=dem, output=fns["dem"], geometry=geom)
        with rasterio.open(fns["dem"]) as f:
            elevation = f.read(1, masked=True)

    # Compute slope
    slope = gdal_slope_analysis(fns["dem"], set_output=fns["slope"])

    # Compute aspect
    aspect = gdal_aspect_analysis(fns["dem"], set_output=fns["aspect"])
    aspect_mean = circular_mean_aspect(aspect)

    return {"elevation": elevation.mean(), "slope": slope.mean(), "aspect": aspect_mean}


def gdal_slope_analysis(
    dem: Union[str, Path],
    set_output: Optional[Union[str, Path]] = None,
    units: str = "degree",
) -> np.ndarray:
    """Return the slope of the terrain from the DEM.

    The slope is the magnitude of the gradient of the elevation.

    Parameters
    ----------
    dem : Union[str, Path]
      Path to file storing DEM.
    set_output : Union[str, Path]
      If set to a valid filepath, will write to this path, otherwise will use an in-memory gdal.Dataset.
    units : str
      Slope units. Default: 'degree'.

    Returns
    -------
    np.ndarray
      Slope array.

    Notes
    -----
    Ensure that the DEM is in a *projected coordinate*, not a geographic coordinate system, so that the
    horizontal scale is the same as the vertical scale (m).

    """
    if isinstance(dem, Path):
        dem = str(dem)
    if set_output:
        if isinstance(set_output, (str, Path)):
            set_output = str(set_output)
            DEMProcessing(
                set_output,
                dem,
                "slope",
                slopeFormat=units,
                format="GTiff",
                band=1,
                creationOptions=[GDAL_TIFF_COMPRESSION_OPTION],
            )
            with rasterio.open(set_output) as src:
                return np.ma.masked_values(src.read(1), value=-9999)
        else:
            raise ValueError()
    else:
        set_output = DEMProcessing(
            "",
            dem,
            "slope",
            slopeFormat=units,
            format="MEM",
            band=1,
        )
        return np.ma.masked_values(set_output.ReadAsArray(), value=-9999)


def gdal_aspect_analysis(
    dem: Union[str, Path],
    set_output: Union[str, Path, bool] = False,
    flat_values_are_zero: bool = False,
) -> Union[np.ndarray, Dataset]:
    """Return the aspect of the terrain from the DEM.

    The aspect is the compass direction of the steepest slope (0: North, 90: East, 180: South, 270: West).

    Parameters
    ----------
    dem : Union[str, Path]
      Path to file storing DEM.
    set_output : Union[str, Path, bool]
      If set to a valid filepath, will write to this path, otherwise will use an in-memory gdal.Dataset.
    flat_values_are_zero: bool
      Designate flat values with value zero. Default: -9999.

    Returns
    -------
    np.ndarray
      Aspect array.

    Notes
    -----
    Ensure that the DEM is in a *projected coordinate*, not a geographic coordinate system, so that the
    horizontal scale is the same as the vertical scale (m).
    """
    if isinstance(dem, Path):
        dem = str(dem)
    if set_output:
        if isinstance(set_output, (str, Path)):
            set_output = str(set_output)
            DEMProcessing(
                destName=set_output,
                srcDS=dem,
                processing="aspect",
                zeroForFlat=flat_values_are_zero,
                format="GTiff",
                band=1,
                creationOptions=[GDAL_TIFF_COMPRESSION_OPTION],
            )
            with rasterio.open(set_output) as src:
                return np.ma.masked_values(src.read(1), value=-9999)
        else:
            raise ValueError()

    else:
        set_output = DEMProcessing(
            destName="",
            srcDS=dem,
            processing="aspect",
            zeroForFlat=flat_values_are_zero,
            format="MEM",
            band=1,
        )
        return np.ma.masked_values(set_output.ReadAsArray(), value=-9999)


def circular_mean_aspect(angles: np.ndarray) -> np.ndarray:
    """Return the mean angular aspect based on circular arithmetic approach

    Parameters
    ----------
    angles: np.ndarray
      Array of aspect angles

    Returns
    -------
    np.ndarray
      Circular mean of aspect array.
    """
    # Circular statistics needed for mean angular aspect
    # Example from: https://gis.stackexchange.com/a/147135/65343

    n = len(angles)
    sine_mean = np.divide(np.sum(np.sin(np.radians(np.ma.masked_array(angles)))), n)
    cosine_mean = np.divide(np.sum(np.cos(np.radians(np.ma.masked_array(angles)))), n)
    vector_mean = np.arctan2(sine_mean, cosine_mean)
    degrees = np.degrees(vector_mean)

    if degrees < 0:
        return degrees + 360
    return degrees


def generic_raster_clip(
    raster: Union[str, Path],
    output: Union[str, Path],
    geometry: Union[Polygon, MultiPolygon, List[Union[Polygon, MultiPolygon]]],
    touches: bool = False,
    fill_with_nodata: bool = True,
    padded: bool = True,
    raster_compression: str = RASTERIO_TIFF_COMPRESSION,
) -> None:
    """
    Crop a raster file to a given geometry.

    Parameters
    ----------
    raster : Union[str, Path]
      Path to input raster.
    output : Union[str, Path]
      Path to output raster.
    geometry : Union[Polygon, MultiPolygon, List[Union[Polygon, MultiPolygon]]
      Geometry defining the region to crop.
    touches : bool
      Whether or not to include cells that intersect the geometry. Default: True.
    fill_with_nodata: bool
      Whether or not to keep pixel values for regions outside of shape or set as nodata. Default: True.
    padded: bool
      Whether or not to add a half-pixel buffer to shape before masking
    raster_compression : str
      Level of data compression. Default: 'lzw'.

    Returns
    -------
    None
    """
    if not (type(geometry) in (list, tuple)):
        geometry = [geometry]

    with rasterio.open(raster, "r") as src:
        mask_image, mask_affine = rasterio.mask.mask(
            src,
            geometry,
            crop=True,
            pad=padded,
            all_touched=touches,
            filled=fill_with_nodata,
        )
        mask_meta = src.meta.copy()
        mask_meta.update(
            {
                "driver": "GTiff",
                "height": mask_image.shape[1],
                "width": mask_image.shape[2],
                "transform": mask_affine,
                "compress": raster_compression,
            }
        )

        # Write the new masked image
        with rasterio.open(output, "w", **mask_meta) as dst:
            dst.write(mask_image)
    return


def generic_raster_warp(
    raster: Union[str, Path],
    output: Union[str, Path],
    target_crs: Union[str, dict, CRS],
    raster_compression: str = RASTERIO_TIFF_COMPRESSION,
) -> None:
    """
    Reproject a raster file.

    Parameters
    ----------
    raster : Union[str, Path]
      Path to input raster.
    output : Union[str, Path]
      Path to output raster.
    target_crs : str or dict
      Target projection identifier.
    raster_compression: str
      Level of data compression. Default: 'lzw'.

    Returns
    -------
    None
    """
    with rasterio.open(raster, "r") as src:
        # Reproject raster using WarpedVRT class
        with rasterio.vrt.WarpedVRT(src, crs=target_crs) as vrt:
            # Calculate grid properties based on projection
            affine, width, height = rasterio.warp.calculate_default_transform(
                src.crs, target_crs, src.width, src.height, *src.bounds
            )

            # Copy relevant metadata from parent raster
            metadata = src.meta.copy()
            metadata.update(
                {
                    "driver": "GTiff",
                    "height": height,
                    "width": width,
                    "transform": affine,
                    "crs": target_crs,
                    "compress": raster_compression,
                }
            )
            data = vrt.read()

            with rasterio.open(output, "w", **metadata) as dst:
                dst.write(data)
    return


def generic_vector_reproject(
    vector: Union[str, Path],
    projected: Union[str, Path],
    source_crs: Union[str, CRS] = WGS84,
    target_crs: Union[str, CRS] = None,
) -> None:
    """Reproject all features and layers within a vector file and return a GeoJSON

    Parameters
    ----------
    vector : Union[str, Path]
      Path to a file containing a valid vector layer.
    projected: Union[str, Path]
      Path to a file to be written.
    source_crs : Union[str, dict, CRS]
      Projection identifier (proj4) for the source geometry, Default: '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, dict, CRS]
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    None
    """

    if target_crs is None:
        raise ValueError("No target CRS is defined.")

    output = {"type": "FeatureCollection", "features": []}

    if isinstance(vector, Path):
        vector = vector.as_posix()

    for i, layer_name in enumerate(fiona.listlayers(vector)):
        with fiona.open(vector, "r", layer=i) as src:
            with open(projected, "w") as sink:
                for feature in src:
                    # Perform vector reprojection using Shapely on each feature
                    try:
                        geom = shape(feature["geometry"])
                        transformed = geom_transform(geom, source_crs, target_crs)
                        feature["geometry"] = mapping(transformed)
                        output["features"].append(feature)
                    except Exception as e:
                        LOGGER.exception(
                            "%s: Unable to reproject feature %s" % (e, feature)
                        )
                        raise

                sink.write(f"{json.dumps(output)}")
    return


def get_bbox(vector: Union[str, Path], all_features: bool = True) -> tuple:
    """Return bounding box of all features or the first feature in file.

    Parameters
    ----------
    vector : str
      Path to file storing vector features.
    all_features : bool
      Return the bounding box for all features. Default: True.

    Returns
    -------
    tuple
      Geographic coordinates of the bounding box (lon0, lat0, lon1, lat1).

    """

    if not all_features:
        with fiona.open(vector, "r") as src:
            for feature in src:
                geom = shape(feature["geometry"])
                return geom.bounds

    with fiona.open(vector, "r") as src:
        return src.bounds


def feature_contains(
    point: Union[Tuple[Union[int, float, str], Union[str, float, int]], Point],
    shp: Union[str, Path, List[Union[str, Path]]],
) -> Union[dict, bool]:
    """Return the first feature containing a location.

    Parameters
    ----------
    point : Union[Tuple[Union[int, float, str], Union[str, float, int]], Point]
      Geographic coordinates of a point (lon, lat) or a shapely Point.
    shp : Union[str, Path, List[str, Path]]
      Path to the file storing the geometries.

    Returns
    -------
    Union[dict, bool]
      The feature found.

    Notes
    -----
    This is really slow. Another approach is to use the `fiona.Collection.filter` method.
    """

    if isinstance(point, collections.abc.Sequence) and not isinstance(point, str):
        for coord in point:
            if isinstance(coord, (int, float)):
                pass
        point = Point(point)
    elif isinstance(point, Point):
        pass
    else:
        raise ValueError(
            f"point should be shapely.Point or tuple of coordinates, got : {point} of type({type(point)})"
        )

    shape_crs = crs_sniffer(single_file_check(shp))

    if isinstance(shp, list):
        shp = shp[0]

    for i, layer_name in enumerate(fiona.listlayers(str(shp))):
        with fiona.open(shp, "r", crs=shape_crs, layer=i) as vector:
            for f in vector.filter(bbox=(point.x, point.y, point.x, point.y)):
                return f

    return False
