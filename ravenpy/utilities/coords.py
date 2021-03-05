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
    model = models.get_model(model)
    return xr.IndexVariable(
        "param",
        data=np.array(model.params._fields),
        attrs={
            "standard_name": "parameter",
            "long_name": "{} model parameter name".format(model),
        },
    )
