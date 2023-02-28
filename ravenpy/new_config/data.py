from .conventions import CF_RAVEN
import xarray as xr


def nc_specs(fn, data_type, station_idx, alt_names=()):
    """Extract specifications from netCDF file.

    Parameters
    ----------
    fn : str, Path
      NetCDF file path.
    data_type: str
      Raven data type.
    station_idx: int
      Index along station dimension. Starts at 1.
    alt_names: list
      Alternative variable names for data type if not the CF standard default.
    """
    from ravenpy.utilities.coords import infer_scale_and_offset

    # Convert to NumPy 0-based indexing
    i = station_idx - 1

    attrs = {"file_name_nc": {"value": fn},
             "data_type": data_type}
    with xr.open_dataset(fn) as ds:
        var_names = [CF_RAVEN[data_type], ] + list(alt_names)
        for v in var_names:
            if v in ds.data_vars:
                nc_var = ds[v]
                attrs["var_name_nc"] = {"value": v}
                attrs["dim_names_nc"] = {"value": nc_var.dims}
                attrs["units"] = nc_var.attrs.get("units")
                if attrs["units"] is not None:
                    s, o = infer_scale_and_offset(nc_var, data_type)
                    attrs["linear_transform"] = dict(scale=s, offset=o)

                break
        else:
            raise ValueError(f"No variable found for {data_type}.\n {ds.data_vars}")

        try:
            attrs["latitude_var_name_nc"] = {"value": ds.cf["latitude"].name}
            attrs["longitude_var_name_nc"] = {"value": ds.cf["longitude"].name}

            attrs["latitude"] = {"value": ds.cf["latitude"][i]}
            attrs["longitude"] = {"value": ds.cf["longitude"][i]}

        except KeyError:
            pass

        try:
            nc_elev = ds.cf["vertical"].name
            attrs["elevation"] = {"value": ds.cf["vertical"][i]}
        except KeyError:
            nc_elev = "elevation" if "elevation" in ds else None
        finally:
            if nc_elev is not None:
                attrs["elevation_var_name_nc"] = {"value": nc_elev}
                attrs["elevation"] = {"value": ds["elevation"][i]}

        if "station_id" in ds:
            if ds["station_id"].shape and len(ds["station_id"]) > i:
                attrs["name"] = {"value": ds["station_id"].values[i]}

        return attrs


def filter_for(kls, attrs):
    """Return attributes that are fields of dataclass."""
    return {k: v for (k,v) in attrs.items() if k in kls.__dataclass_fields__.keys()}
