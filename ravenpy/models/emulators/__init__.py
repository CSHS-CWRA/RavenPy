from .blended import *
from .canadianshield import *
from .gr4jcn import *
from .hbvec import *
from .hmets import *
from .hypr import *
from .mohyse import *
from .sacsma import *


def get_model(name):
    """Return the corresponding Raven emulated model instance.

    Parameters
    ----------
    name : str
      Model class name or model identifier.

    Returns
    -------
    Raven model instance
    """
    from ravenpy.models import emulators

    model_cls = getattr(emulators, name.upper(), None)

    if model_cls is None:
        raise ValueError("Model {} is not recognized.".format(name))

    return model_cls
