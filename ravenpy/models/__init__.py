import os

from .aggregate_forcings_to_HRUs import aggregate_forcings_to_HRUs
from .base import Ostrich, Raven
from .emulators import *
from .multimodel import RavenMultiModel
from .rv import RV, RVI

_dir = os.path.abspath(os.path.dirname(__file__))

raven_templates = {
    "raven-gr4j-cemaneige": os.path.join(_dir, "raven-gr4j-cemaneige"),
    "raven-mohyse": os.path.join(_dir, "raven-mohyse"),
    "raven-hmets": os.path.join(_dir, "raven-hmets"),
    "raven-hbv-ec": os.path.join(_dir, "raven-hbv-ec"),
}

ostrich_templates = {
    "ostrich-gr4j-cemaneige": os.path.join(_dir, "ostrich-gr4j-cemaneige"),
    "ostrich-mohyse": os.path.join(_dir, "ostrich-mohyse"),
    "ostrich-hmets": os.path.join(_dir, "ostrich-hmets"),
    "ostrich-hbv-ec": os.path.join(_dir, "ostrich-hbv-ec"),
}
