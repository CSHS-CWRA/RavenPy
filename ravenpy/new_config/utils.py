import dataclasses

import xarray as xr

from .conventions import CF_RAVEN, MonthlyAverages


def nc_specs(fn, data_type, station_idx, alt_names=(), mon_ave=False):
    """Extract specifications from netCDF file.

    Parameters
    ----------
    fn : str, Path
      NetCDF file path.
    data_type: str
      Raven data type.
    station_idx: int
      Index along station dimension. Starts at 1.
    alt_names: str, list
      Alternative variable names for data type if not the CF standard default.
    mon_ave: bool
      If True, compute the monthly average.

    Returns
    -------
    dict

    Notes
    -----
    Attributes: var_name_nc, dim_names_nc, units, linear_transform, latitude_var_name_nc, longitude_var_name_nc,
    elevation_var_name_nc
    latitude, longitude, elevation, name

    """
    from pathlib import Path

    from ravenpy.utilities.coords import infer_scale_and_offset

    # Convert to NumPy 0-based indexing
    i = station_idx - 1

    if Path(fn).exists():
        fn = Path(fn).resolve(strict=True)
    else:
        raise ValueError("NetCDF file not found.")

    if isinstance(alt_names, str):
        alt_names = (alt_names,)

    attrs = {"file_name_nc": fn, "data_type": data_type, "forcing_type": data_type}

    with xr.open_dataset(fn) as ds:
        var_names = (CF_RAVEN.get(data_type),) + alt_names
        for v in var_names:
            if v in ds.data_vars:
                nc_var = ds[v]
                attrs["var_name_nc"] = v
                attrs["dim_names_nc"] = nc_var.dims
                attrs["_time_dim_name_nc"] = ds.cf["time"].name
                attrs["_dim_size_nc"] = dict(zip(nc_var.dims, nc_var.shape))
                attrs["units"] = nc_var.attrs.get("units")
                if attrs["units"] is not None:
                    s, o = infer_scale_and_offset(nc_var, data_type)
                    attrs["linear_transform"] = dict(scale=s, offset=o)
                if mon_ave:
                    ma = MonthlyAverages.get(data_type)
                    if ma:
                        attrs[ma] = nc_var.groupby("time.month").mean().values

                break
        else:
            raise ValueError(f"No variable found for {data_type}.\n {ds.data_vars}")

        try:
            attrs["latitude_var_name_nc"] = ds.cf["latitude"].name
            attrs["longitude_var_name_nc"] = ds.cf["longitude"].name

            attrs["latitude"] = ds.cf["latitude"][i]
            attrs["longitude"] = ds.cf["longitude"][i]

        except KeyError:
            pass

        try:
            nc_elev = ds.cf["vertical"].name
            attrs["elevation"] = ds.cf["vertical"][i]
        except KeyError:
            nc_elev = "elevation" if "elevation" in ds else None
        finally:
            if nc_elev is not None:
                attrs["elevation_var_name_nc"] = nc_elev
                attrs["elevation"] = ds["elevation"][i]

        if "station_id" in ds:
            if ds["station_id"].shape and len(ds["station_id"]) > i:
                attrs["name"] = ds["station_id"].values[i]

        return attrs


def filter_for(kls, attrs):
    """Return attributes that are fields of dataclass."""
    from pydantic.main import ModelMetaclass

    from .commands import Command

    if hasattr(kls, "__dataclass_fields__"):
        return {
            k: v for (k, v) in attrs.items() if k in kls.__dataclass_fields__.keys()
        }
    else:
        out = {}
        for key, field in kls.__fields__.items():
            if key in attrs:
                out[key] = attrs[key]
            elif type(field.type_) == ModelMetaclass and issubclass(
                field.type_, Command
            ):
                out[key] = filter_for(field.type_, attrs)
        return out
        # return {k: v for (k, v) in attrs.items() if k in kls.__fields__.keys()}
