from dataclasses import fields

import numpy as np
import xarray as xr

import ravenpy.models as models


def realization(n):
    """Return a realization coordinate.

    Parameters
    ----------
    n : int
      Size of the ensemble.
    """
    return xr.IndexVariable(
        "realization",
        data=range(n),
        attrs={
            "standard_name": "realization",
            "axis": "E",
            "units": "1",
            "long_name": "Label identifying the ensemble member",
        },
    )


def param(model):
    """Return a parameter coordinate.

    Parameters
    ----------
    model : str
      Model name.
    """
    model_cls = models.get_model(model)
    return xr.IndexVariable(
        "param",
        data=np.array([f.name for f in fields(model_cls.Params)]),
        attrs={
            "standard_name": "parameter",
            "long_name": f"{model} model parameter name",
        },
    )


def infer_scale_and_offset(da: xr.DataArray, data_type: str) -> (float, float):
    """Return scale and offset parameters from data.

    Infer scale and offset parameters describing the linear transformation from the units in file to Raven
    compliant units.

    Parameters
    ----------
    da : xr.DataArray
      Input data.
    data_type : str
      Raven data type, e.g. 'PRECIP', 'TEMP_AVE', etc.

    Returns
    -------
    float, float
      Scale and offset parameters.
    """
    import pint
    from xclim.core.units import FREQ_UNITS, parse_offset, units, units2pint

    from ravenpy.config import defaults

    # Read units attribute from netCDF file -- assuming CF-Compliance.
    source = units2pint(da.attrs["units"])

    # Get default units for data type.
    target = defaults.units[data_type]

    # Linear transform parameters
    try:
        scale, offset = units_transform(source, target)
    except pint.errors.DimensionalityError:
        if data_type in ["PRECIP", "PRECIP_DAILY_AVE", "RAINFALL", "SNOWFALL"]:
            # Source units are in total precipitation instead of rate. We need to infer the accumulation time
            # in order to find the transform.
            freq = xr.infer_freq(da.time)
            if freq is None:
                raise ValueError(f"Cannot infer time frequency of input data {da}")
            multi, base, start_anchor, _ = parse_offset(freq)
            if base in ["M", "Q", "A"]:
                raise ValueError(f"Irregular time frequency for input data {da}")
            source = source / multi / units(FREQ_UNITS[base])
            scale, offset = units_transform(source, target)

    return scale, offset


def units_transform(source, target):
    """Return linear transform parameters to convert one unit to another.

    If the target unit is given by `y = ax + b`, where `x` is the value of the source unit, then this function
    returns a, b.

    Parameters
    ----------
    source : str, pint.Unit
        Source unit string, pint-recognized.
    target : str
        Target unit string, pint-recognized.
    """
    from xclim.core.units import convert_units_to, units

    b = convert_units_to(units.Quantity(0, source), target)
    a = convert_units_to(units.Quantity(1, source), target) - b

    return a, b
