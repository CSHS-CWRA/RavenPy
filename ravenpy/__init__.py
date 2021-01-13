"""Top-level package for RavenPy."""

import os
import shutil

from .__version__ import __author__, __email__, __version__  # noqa: F401

RAVEN_EXEC_PATH = os.getenv("RAVENPY_RAVEN_BINARY_PATH") or shutil.which("raven")
OSTRICH_EXEC_PATH = os.getenv("RAVENPY_OSTRICH_BINARY_PATH") or shutil.which("ostrich")

if not RAVEN_EXEC_PATH:
    raise RuntimeError(
        "Could not find raven binary in PATH, and RAVEN_BINARY_PATH env variable is not set"
    )

if not OSTRICH_EXEC_PATH:
    raise RuntimeError(
        "Could not find ostrich binary in PATH, and OSTRICH_BINARY_PATH env variable is not set"
    )
