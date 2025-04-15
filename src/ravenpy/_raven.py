"""Configurations for raven-hydro."""

import os
import shutil

RAVEN_EXEC_PATH = os.getenv("RAVENPY_RAVEN_BINARY_PATH") or shutil.which("raven")

if not RAVEN_EXEC_PATH:
    raise RuntimeError(
        "Could not find raven binary in PATH and RAVENPY_RAVEN_BINARY_PATH env variable is not set."
    )

try:
    from raven_hydro import __raven_version__
except ImportError:
    __raven_version__ = "0.0.0"
