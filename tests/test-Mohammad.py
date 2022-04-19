import warnings
from collections import defaultdict
from pathlib import Path
from typing import List

from ravenpy.utilities import gis_import_error_message

try:
    import geopandas
    from osgeo import __version__ as osgeo_version  # noqa
    from osgeo import ogr, osr  # noqa
    from shapely import wkt
except (ImportError, ModuleNotFoundError) as e:
    msg = gis_import_error_message.format(Path(__file__).stem)
    raise ImportError(msg) from e

import netCDF4 as nc4
import numpy as np

from ravenpy.config.commands import (
    ChannelProfileCommand,
    GridWeightsCommand,
    HRUsCommand,
    ReservoirCommand,
    SubBasinsCommand,
)
from ravenpy.extractors.routing_product import RoutingProductShapefileExtractor
from ravenpy.utilities.testdata import get_local_testdata

# routing_product_shp_path = get_local_testdata("Famine/HRU_Famine.zip")
pth = "/home/mohammad/Dossier_travail/Hydrotel/DEH/MG24HA/SLSO_MG24HA_2020/physitel/HRU/raven-testdata/famine/hru_Famine_final.zip"

rvh_extractor = RoutingProductShapefileExtractor(
    pth,
    hru_aspect_convention="ArcGIS",
    routing_product_version="2.1",
)
rvh_config = rvh_extractor.extract()

###################################################################################################################################

import collections
import csv
import datetime as dt
import operator
import os
import re
import shutil
import stat
import subprocess
import tempfile
import zipfile
from collections import OrderedDict
from dataclasses import astuple, fields, is_dataclass, replace
from pathlib import Path
from typing import Any, Dict, List, Union, cast
from warnings import warn

import numpy as np
import xarray as xr
from numpy.distutils.misc_util import is_sequence

from ravenpy.extractors.routing_product import RoutingProductShapefileExtractor
from ravenpy.models import Raven

model = Raven(identifier="famine", workdir="./test_rv_Mohammad")
