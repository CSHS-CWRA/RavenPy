"""A Python wrapper for configuring and running the hydrologic modelling framework Raven."""

###################################################################################
# MIT License
#
# Copyright (c) 2020-2025 David Huard, Trevor James Smith, Christian Jauvin, Julie Mai, Ming Han
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###################################################################################

import os
import shutil

from .ravenpy import Emulator, EnsembleReader, OutputReader, RavenWarning, run

__all__ = [
    "RAVEN_EXEC_PATH",
    "Emulator",
    "EnsembleReader",
    "OutputReader",
    "RavenWarning",
    "__author__",
    "__email__",
    "__raven_version__",
    "__version__",
    "run",
]


RAVEN_EXEC_PATH = os.getenv("RAVENPY_RAVEN_BINARY_PATH") or shutil.which("raven")

if not RAVEN_EXEC_PATH:
    raise RuntimeError(
        "Could not find raven binary in PATH and RAVENPY_RAVEN_BINARY_PATH env variable is not set."
    )

try:
    from raven_hydro import __raven_version__
except ImportError:
    __raven_version__ = "0.0.0"

__author__ = """David Huard"""
__email__ = "huard.david@ouranos.ca"
__version__ = "0.18.0"
