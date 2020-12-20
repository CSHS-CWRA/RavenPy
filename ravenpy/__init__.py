"""Top-level package for RavenPy."""

import os
import sys
import warnings
from pathlib import Path

from .__version__ import __author__, __email__, __version__  # noqa: F401

# Note: all the binaries associated with this package are installed in the "bin"
#       folder of the target venv, which corresponds to either `sys.prefix` or
#       $CONDA_PREFIX.


if os.getenv("CONDA_PREFIX"):
    VENV_PATH = Path(os.getenv("CONDA_PREFIX"))
elif sys.base_prefix != sys.prefix:
    VENV_PATH = Path(sys.prefix)
else:
    raise RuntimeError("Please run RavenPy in a virtual environment!")


RAVEN_EXEC_PATH = VENV_PATH / "bin" / "raven"
OSTRICH_EXEC_PATH = VENV_PATH / "bin" / "ostrich"


# if "DO_NOT_CHECK_EXECUTABLE_EXISTENCE" not in os.environ:
#     raven_exec = VENV_PATH / "bin" / "raven"
#     if not raven_exec.exists():
#         raise IOError(
#             "The raven executable is not installed in the virtual environment."
#         )

#     ostrich_exec = VENV_PATH / "bin" / "ostrich"
#     if not ostrich_exec.exists():
#         raise IOError(
#             "The ostrich executable is not installed in the virtual environment."
#         )

# raven_simg = VENV_PATH / "bin" / "hydro-raven-latest.simg"
# if not raven_simg.exists():
#     warnings.warn(
#         "The Raven Singularity image has not been downloaded. Execute \n"
#         "$ singularity pull shub://132.217.141.54/hydro/raven:latest \n"
#         "and store the image in <your_venv>/bin/"
#     )
