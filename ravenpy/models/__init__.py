import os

from .base import Ostrich, Raven
from .emulators import *
from .multimodel import RavenMultiModel
from .rv import HRU, LU, RV, RVI, HRUState

_dir = os.path.abspath(os.path.dirname(__file__))

raven_templates = {
    "raven-gr4j-cemaneige": os.path.join(_dir, "raven-gr4j-cemaneige"),
    "raven-mohyse": os.path.join(_dir, "raven-mohyse"),
    "raven-hmets": os.path.join(_dir, "raven-hmets"),
    "raven-hbv-ec": os.path.join(_dir, "raven-hbv-ec"),
    "raven-blended": os.path.join(_dir, "raven-blended"),
}

ostrich_templates = {
    "ostrich-gr4j-cemaneige": os.path.join(_dir, "ostrich-gr4j-cemaneige"),
    "ostrich-mohyse": os.path.join(_dir, "ostrich-mohyse"),
    "ostrich-hmets": os.path.join(_dir, "ostrich-hmets"),
    "ostrich-hbv-ec": os.path.join(_dir, "ostrich-hbv-ec"),
    "ostrich-blended": os.path.join(_dir, "ostrich-blended"),
}

grid_weight_importer_params = dict(
    DIM_NAMES=("lon_dim", "lat_dim"),
    VAR_NAMES=("lon", "lat"),
    ROUTING_ID_FIELD="HRU_ID",
    NETCDF_INPUT_FIELD="NetCDF_col",
    AREA_ERROR_THRESHOLD=0.05,
)
