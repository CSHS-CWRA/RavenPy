from ravenpy.config.emulators.blended import Blended
from ravenpy.config.emulators.canadianshield import CanadianShield
from ravenpy.config.emulators.gr4jcn import GR4JCN
from ravenpy.config.emulators.hbvec import HBVEC
from ravenpy.config.emulators.hmets import HMETS
from ravenpy.config.emulators.hypr import HYPR
from ravenpy.config.emulators.mohyse import Mohyse
from ravenpy.config.emulators.routing import BasicRoute
from ravenpy.config.emulators.sacsma import SACSMA

__all__ = [
    "get_model",
]


def get_model(name):
    """
    Return the corresponding Raven emulator configuration class.

    Parameters
    ----------
    name : str
        Model class name or model identifier.

    Returns
    -------
    Raven model configuration class
    """
    cls = globals().get(name)

    if cls is None:
        raise ValueError(f"Model {name} is not recognized.")

    return cls
