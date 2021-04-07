from .blended import *
from .gr4jcn import *
from .hbvec import *
from .hmets import *
from .mohyse import *


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

    model_cls = getattr(emulators, name, None)

    if model_cls is None:
        for m in [GR4JCN, MOHYSE, HMETS, HBVEC, BLENDED]:
            if m.identifier == name:
                model_cls = m

    if model_cls is None:
        raise ValueError("Model {} is not recognized.".format(name))

    return model_cls
