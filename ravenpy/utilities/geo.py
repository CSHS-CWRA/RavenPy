"""
Tools for performing geospatial translations and transformations.
"""

import collections
import json
import logging
from pathlib import Path
from typing import List, Union

from . import gis_import_error_message

try:
    import fiona
    import rasterio
    import rasterio.mask
    import rasterio.vrt
    import rasterio.warp
    from pyproj import CRS
    from shapely.geometry import (
        GeometryCollection,
        MultiPolygon,
        Polygon,
        mapping,
        shape,
    )
    from shapely.ops import transform
except (ImportError, ModuleNotFoundError) as e:
    msg = gis_import_error_message.format(Path(__file__).stem)
    raise ImportError(msg) from e

RASTERIO_TIFF_COMPRESSION = "lzw"
LOGGER = logging.getLogger("RavenPy")
WGS84 = 4326


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
        from functools import partial

        from pyproj import Transformer  # noqa

        if isinstance(source_crs, int or str):
            source = CRS.from_epsg(source_crs)
        else:
            source = source_crs

        if isinstance(target_crs, int or str):
            target = CRS.from_epsg(target_crs)
        else:
            target = target_crs

        transform_func = Transformer.from_crs(source, target, always_xy=True)
        reprojected = transform(transform_func.transform, geom)

        return reprojected
    except Exception as err:
        msg = f"{err}: Failed to reproject geometry"
        LOGGER.error(msg)
        raise Exception(msg)


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
    if not isinstance(geometry, collections.abc.Iterable):
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

    output = {"type": "FeatureCollection", "features": list()}

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
                    except Exception as err:
                        LOGGER.exception(
                            "{}: Unable to reproject feature {}".format(err, feature)
                        )
                        raise

                sink.write(f"{json.dumps(output)}")
