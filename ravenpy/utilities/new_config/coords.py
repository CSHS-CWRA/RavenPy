from dataclasses import fields

import numpy as np
import xarray as xr

from ravenpy.new_config.emulators import get_model


def param(model):
    """Return a parameter coordinate.

    Parameters
    ----------
    model : str
      Model name.
    """
    cls = get_model(model)
    P = cls.__fields__["params"].type_
    return xr.IndexVariable(
        "param",
        data=np.array([f.name for f in fields(P)]),
        attrs={
            "standard_name": "parameter",
            "long_name": f"{model} model parameter name",
        },
    )
