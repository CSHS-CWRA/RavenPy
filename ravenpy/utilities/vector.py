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
from os import PathLike
from pathlib import Path
from typing import Iterable, List, Optional, Union

import geojson
import pandas as pd
import pygml
import shapefile
from pyproj import CRS
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


def geom_properties(geom: Union[Polygon, MultiPolygon, GeometryCollection]) -> dict:
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
) -> GeometryCollection:
    """Change the projection of a geometry.

    Assuming a geometry's coordinates are in a `source_crs`, compute the new coordinates under the `target_crs`.

    Parameters
    ----------
    geom : shapely.geometry.GeometryCollection or shapely.geometry.shape
      Source geometry.
    source_crs : int or str or pyproj.CRS
      Projection identifier (proj4) for the source geometry, e.g. '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : int or str or pyproj.CRS
      Projection identifier (proj4) for the target geometry.
    always_xy : bool

    Returns
    -------
    shapely.geometry.GeometryCollection
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

        transform_func = Transformer.from_crs(source, target, always_xy=always_xy)
        reprojected = transform(transform_func.transform, geom)

        return reprojected
    except Exception as err:
        msg = f"{err}: Failed to reproject geometry."
        LOGGER.error(msg)
        raise Exception(msg)


def geojson_object_transform(
    collection: geojson.FeatureCollection,
    source_crs: Union[int, str, CRS],
    target_crs: Union[int, str, CRS],
    always_xy: bool = True,
) -> geojson.FeatureCollection:
    """

    Parameters
    ----------
    collection : geojson.FeatureCollection
    source_crs : int or str or pyproj.CRS
    target_crs : int or str or pyproj.CRS
    always_xy : bool

    Returns
    -------
    geojson.FeatureCollection
    """
    output = list()
    for feature in collection.features:
        try:
            geom = shape(feature.geometry)
            transformed = geom_transform(
                geom, source_crs, target_crs, always_xy=always_xy
            )
            feature.geometry = mapping(transformed)
            if hasattr(feature, "bbox"):
                bbox = box(*feature.bbox)
                transformed = geom_transform(
                    bbox, source_crs, target_crs, always_xy=always_xy
                )
                feature.bbox = mapping(transformed)
            output.append(feature)

        except Exception as err:
            LOGGER.exception(
                "{}: Unable to reproject feature {}".format(err, collection)
            )
            raise

    return geojson.FeatureCollection(output)


def generic_vector_file_transform(
    vector: Union[str, PathLike],
    projected: Union[str, PathLike],
    source_crs: Union[str, CRS] = WGS84,
    target_crs: Union[str, CRS] = None,
) -> Union[str, PathLike]:
    """Reproject all features and layers within a vector file and return a GeoJSON

    Parameters
    ----------
    vector : str or PathLike
      Path to a file containing a valid vector layer (GeoJSON, Shapefile, or vector GML).
    projected: str or PathLike
      Path to a file to be written.
    source_crs : Union[str, dict, CRS]
      Projection identifier (proj4) for the source geometry, Default: '+proj=longlat +datum=WGS84 +no_defs'.
    target_crs : Union[str, dict, CRS]
      Projection identifier (proj4) for the target geometry.

    Returns
    -------
    str or PathLike
    """

    if target_crs is None:
        raise ValueError("No target CRS is defined.")

    if Path(vector).suffix in [".tar", ".zip", ".7z"]:
        warnings.warn(
            f"File {Path(vector).name} should be extracted from archive before reprojection."
        )
        vector = next(iter(archive_sniffer(_generic_extract_archive(vector))))

    if isinstance(vector, Path):
        vector = vector.as_posix()

    if vector.lower().endswith(".shp"):
        src = geojson.loads(json.dumps(shapefile.Reader(vector).__geo_interface__))
    elif vector.lower().endswith("json"):
        src = geojson.load(open(vector))
    elif vector.lower().endswith("gml"):
        src = geojson.load(
            json.dumps(pygml.parse(open(vector, "rb").read()).__geo_interface__)
        )
    else:
        raise FileNotFoundError(f"{vector} is not a valid GeoJSON or Shapefile.")

    with open(projected, "w") as sink:
        output = geojson_object_transform(src, source_crs, target_crs)
        sink.write(f"{json.dumps(output)}")

    return projected


def _generic_extract_archive(
    resources: Union[str, PathLike, List[Union[bytes, str, PathLike]]],
    output_dir: Optional[Union[str, PathLike]] = None,
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
    archives: Union[str, PathLike, Iterable[Union[str, PathLike]]],
    working_dir: Optional[Union[str, PathLike]] = None,
    extensions: Optional[Iterable[str]] = None,
) -> List[Union[str, PathLike]]:
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

    decompressed_files = _generic_extract_archive(archives, output_dir=working_dir)
    for file in decompressed_files:
        if any(ext in Path(file).suffix for ext in extensions):
            potential_files.append(file)
    return potential_files


def vector_to_dataframe(
    shp_path: Union[str, PathLike, shapefile.Shape, geojson.feature.FeatureCollection]
) -> pd.DataFrame:
    """
    Read a shapefile into a Pandas dataframe with a 'coords' column holding
    the geometry information. This uses the pyshp package.

    Parameters
    ----------
    shp_path : str or os.PathLike or shapefile.Shape or geojson.feature.FeatureCollection

    Notes
    -----
    Example adapted from https://gist.github.com/aerispaha/f098916ac041c286ae92d037ba5c37ba
    """
    if isinstance(shp_path, (str, PathLike)):
        shp_path = Path(shp_path)

    try:
        import geopandas as gpd

        df = gpd.read_file(shp_path)

    except ModuleNotFoundError:
        try:
            if isinstance(shp_path, shapefile.Shape):
                sf = shp_path
            else:
                sf = shapefile.Reader(shp_path.as_posix())

            fields = [x[0] for x in sf.fields][1:]
            records = sf.records()
            geoms = [s for s in sf.shapes()]

            # write into a dataframe
            df = pd.DataFrame(columns=fields, data=records)
            df = df.assign(geometry=geoms)
            df.geometry = df.geometry.apply(shape)

        except shapefile.ShapefileException:
            try:
                if shp_path.suffix.lower() in {"json", "geojson"}:
                    shp_path = geojson.load(open(shp_path))

                if isinstance(shp_path, geojson.feature.FeatureCollection):
                    df = pd.json_normalize(shp_path.features)
                    df.geometry = df.df.geometry.apply(shape)
                else:
                    raise FileNotFoundError()

            except (FileNotFoundError, UnicodeDecodeError) as e:
                raise ValueError(
                    "Unable to find vector in {}".format(str(shp_path))
                ) from e

    return df
