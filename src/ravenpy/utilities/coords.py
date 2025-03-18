from dataclasses import fields
from typing import Union

import numpy as np
import pint
import xarray as xr

from ravenpy.config.emulators import get_model


def realization(n):
    """Return a realization coordinate.

    Parameters
    ----------
    n : int
        Size of the ensemble.
    """
    return xr.IndexVariable(
        "members",
        data=range(n),
        attrs={
            "standard_name": "members",
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
    cls = get_model(model)
    P = cls.model_fields["params"].annotation  # noqa: N806
    return xr.IndexVariable(
        "param",
        data=np.array([f.name for f in fields(P)]),
        attrs={
            "standard_name": "parameter",
            "long_name": f"{model} model parameter name",
        },
    )


# FIXME: cumulative is not used
def infer_scale_and_offset(
    da: xr.DataArray, data_type: str, cumulative: bool = False  # noqa: F841
) -> tuple[float, float]:
    """
    Return scale and offset parameters from data.

    Infer scale and offset parameters describing the linear transformation from the units in file to Raven
    compliant units.

    Parameters
    ----------
    da : xr.DataArray
        Input data.
    data_type : str
        Raven data type, e.g. 'PRECIP', 'TEMP_AVE', etc.
    cumulative : bool
        Default: False.

    Returns
    -------
    float, float
        Scale and offset parameters.

    Notes
    -----
    Does not work with accumulated variables.
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
            real_source = source / multi / units(FREQ_UNITS.get(base, base))
            scale, offset = units_transform(real_source, target)
        else:
            raise NotImplementedError(f"data_type: {data_type}")

    return scale, offset


def units_transform(source: Union[str, pint.Unit], target: str, context: str = "hydro"):
    """Return linear transform parameters to convert one unit to another.

    If the target unit is given by `y = ax + b`, where `x` is the value of the source unit, then this function
    returns a, b.

    Parameters
    ----------
    source : str, pint.Unit
        Source unit string, pint-recognized.
    target : str
        Target unit string, pint-recognized.
    context : str
        Context of unit conversion. Default: "hydro".
    """
    from xclim.core.units import convert_units_to

    b = convert_units_to(0 * source, target, context)
    a = convert_units_to(1 * source, target, context) - b

    return round(a, 9), round(b, 9)
