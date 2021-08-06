"""
Tools for performing geospatial translations and transformations.
"""
import json
import logging
import math
import tarfile
import tempfile
import warnings
import zipfile
from pathlib import Path
from typing import Iterable, List, Optional, Union

import geojson
from pyproj import CRS
from shapefile import Reader
from shapely.geometry import (
    GeometryCollection,
    MultiPolygon,
    Polygon,
    box,
    mapping,
    shape,
)
from shapely.ops import transform

RASTERIO_TIFF_COMPRESSION = "lzw"
LOGGER = logging.getLogger("RavenPy")
WGS84 = 4326


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


def geom_transform(
    geom: Union[GeometryCollection, shape],
    source_crs: Union[str, int, CRS] = WGS84,
    target_crs: Union[str, int, CRS] = None,
    always_xy: bool = True,
    **kwargs,
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
    always_xy : bool

    **kwargs
      Keyword arguments passed directly to pyproj.Transformer().

    Returns
    -------
    GeometryCollection
      Reprojected geometry.
    """
    try:
        from functools import partial

        from pyproj import Transformer  # noqa

        if isinstance(source_crs, int or str):
            source = CRS.from_user_input(source_crs)
        else:
            source = source_crs

        if isinstance(target_crs, int or str):
            target = CRS.from_user_input(target_crs)
        else:
            target = target_crs

        transform_func = Transformer.from_crs(
            source, target, always_xy=always_xy, **kwargs
        )
        reprojected = transform(transform_func.transform, geom)

        return reprojected
    except Exception as err:
        msg = f"{err}: Failed to reproject geometry."
        LOGGER.error(msg)
        raise Exception(msg)


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
      Path to a file containing a valid vector layer (geoJSON or Shapefile).
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

    output = {"type": "FeatureCollection", "features": list()}

    if Path(vector).suffix in [".tar", ".zip", ".7z"]:
        warnings.warn(
            f"File {Path(vector).name} should be extracted from archive before reprojection."
        )
        vector = next(iter(archive_sniffer(generic_extract_archive(vector))))

    if isinstance(vector, Path):
        vector = vector.as_posix()

    if vector.lower().endswith(".shp"):
        src = geojson.loads(json.dumps(Reader(vector).__geo_interface__))
    elif vector.lower().endswith("json"):
        src = geojson.load(open(vector))
    else:
        raise FileNotFoundError(f"{vector} is not a valid geoJSON or Shapefile.")

    with open(projected, "w") as sink:
        for feature in src.features:
            # Perform vector reprojection using Shapely on each feature
            try:
                geom = shape(feature.geometry)
                transformed = geom_transform(geom, source_crs, target_crs)
                feature.geometry = mapping(transformed)
                if hasattr(feature, "bbox"):
                    bbox = box(*feature.bbox)
                    transformed = geom_transform(bbox, source_crs, target_crs)
                    feature.bbox = mapping(transformed)
                output["features"].append(feature)
            except Exception as err:
                LOGGER.exception("{}: Unable to reproject feature {}".format(err, src))
                raise

        sink.write(f"{json.dumps(output)}")


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
                LOGGER.error(
                    "Failed to extract sub archive {{{}}}: {{{}}}".format(arch, e)
                )
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
