"""Main module."""

import collections
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Union
from warnings import warn

from ravenpy.new_config.rvs import RVC, Config

RAVEN_EXEC_PATH = os.getenv("RAVENPY_RAVEN_BINARY_PATH") or shutil.which("raven")


class Emulator:
    def __init__(self, config: Config, path: Union[Path, str]):
        """
        Convenience class to work with the Raven modeling framework.

        Parameters
        ----------
        config : Config
          Emulator Config instance fully parameterized, i.e. without symbolic expressions.
        path : str, Path
          Path to rv files and model outputs.
        """
        self._config = config.copy()
        self._path = Path(path)

        # Emulator config is made immutable to avoid complex behavior.
        self._config.__config__.allow_mutation = False

    def write_rv(self, overwrite=False):
        """Write the configuration files to disk."""
        self._config.write_rv(workdir=self.path, overwrite=overwrite)

    def run(self):
        """Run the model."""
        self._outputdir = run(self.config.run_name, self.path)
        return self._outputdir

    def parse_outputs(self):
        """Return model outputs."""
        return parse_outputs(self.config.run_name, self._outputdir)

    def parse_diagnostics(self):
        """Return model diagnostics."""
        return parse_diagnostics(
            self._outputdir / f"{self.config.run_name}_Diagnostics.csv"
        )

    @property
    def config(self):
        """Read-only model configuration."""
        return self._config

    @property
    def path(self):
        """Path to RV files and output sub-directory."""
        return self._path


def run(
    run_name: str,
    configdir: Union[str, Path],
    outputdir: Union[str, Path] = None,
    overwrite: bool = True,
):
    """
    Run Raven given the path to a model configuration.

    Parameters
    ----------
    run_name : str
      Configuration files stem, i.e. the file name without extension.
    configdir : str, Path
      Path to configuration files directory.
    outputdir: str, Path
      Path to model simulation output. If None, will write to configdir/output.

    Return
    ------
    Path
      Path to model outputs.

    """
    if not RAVEN_EXEC_PATH:
        raise RuntimeError(
            "Could not find raven binary in PATH, and RAVENPY_RAVEN_BINARY_PATH env variable is not set"
        )

    # Confirm configdir exists
    configdir = Path(configdir)
    if not configdir.exists():
        raise OSError("Workdir should include configuration files.")

    # Create outputdir
    outputdir = Path(outputdir or "output")
    if not outputdir.is_absolute():
        outputdir = configdir / outputdir

    if overwrite and outputdir.exists():
        shutil.rmtree(str(outputdir))

    if not outputdir.exists():
        os.makedirs(str(outputdir))

    # Launch executable, wait for completion.
    cmd = [RAVEN_EXEC_PATH, run_name, "-o", str(outputdir)]

    process = subprocess.Popen(
        cmd,
        cwd=configdir,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    )

    stdout, stderr = process.communicate(input="\n")
    returncode = process.wait()

    if returncode != 0:
        raise OSError(f"Raven segfaulted : \n{stdout}")

    print(process.returncode)
    # Deal with errors and warnings
    messages = extract_raven_messages(configdir)

    for msg in messages["WARNING"] + messages["ADVISORY"]:
        warn(msg, category=RavenWarning)

    if messages["ERROR"]:
        raise RavenError("\n".join(messages["ERROR"]))

    return outputdir


def parse_outputs(run_name, outputdir: [str, Path]):
    """Parse outputs from model execution.

    Parameters
    ----------
    run_name: str
      RunName value identifying model outputs.
    outputdir: str, Path
      Path to model output directory.

    Returns
    -------
    dict
      Dictionary holding model outputs:
        - hydrograph: xarray.Dataset
        - storage: xarray.Dataset
        - solution: RVC
        - diagnostics: Dict[str, list]
      Values are set to None if no file is found.

    """
    if type(outputdir) == str:
        outputdir = Path(outputdir)

    out = {}
    out["solution"] = parse_solution(outputdir / f"{run_name}_solution.rvc")
    out["hydrograph"] = parse_nc(outputdir / f"{run_name}_Hydrographs.nc")
    out["storage"] = parse_nc(outputdir / f"{run_name}_WatershedStorage.nc")
    out["diagnostics"] = parse_diagnostics(outputdir / f"{run_name}_Diagnostics.csv")

    return out


def parse_diagnostics(fn: Path):
    """Return dictionary of performance metrics."""
    import csv

    if fn.exists():
        out = collections.defaultdict(list)

        with open(fn) as f:
            reader = csv.reader(f.readlines())
            header = next(reader)
            for row in reader:
                for key, val in zip(header, row):
                    if "DIAG" in key:
                        val = float(val)  # type: ignore
                    out[key].append(val)

            out.pop("")

        return out


def parse_nc(fn: Path):
    """Open netCDF dataset with xarray if the path is valid, otherwise return None."""
    import xarray as xr

    if fn.exists():
        if fn.suffix == ".nc":
            return xr.open_dataset(fn)
        else:
            raise NotImplementedError


def parse_solution(fn: Path):
    """Create RVC from model output."""
    if fn.exists():
        if fn.suffix == ".rvc":
            return RVC.from_solution(fn.read_text())
        else:
            raise NotImplementedError


def extract_raven_messages(path):
    """Parse Raven_errors and extract the messages, structured by types."""
    err_filepaths = path.rglob("Raven_errors.txt")
    messages = {
        "ERROR": [],
        "WARNING": [],
        "ADVISORY": [],
        "SIMULATION COMPLETE": False,
    }
    for p in err_filepaths:
        # The error message for an unknown command is exceptionally on two lines
        # (the second starts with a triple space)
        for m in re.findall("^([A-Z ]+) :(.+)(?:\n   (.+))?", p.read_text(), re.M):
            if m[0] == "SIMULATION COMPLETE":
                messages["SIMULATION COMPLETE"] = True
                continue
            msg_type = m[0]
            msg = f"{m[1]} {m[2]}".strip()
            if msg == "Errors found in input data. See Raven_errors.txt for details":
                # Skip this one because it's a bit circular
                continue
            messages[msg_type].append(msg)  # type: ignore

    return messages


class RavenError(Exception):
    """
    This is an error that is meant to be raised whenever a message of type "ERROR" is found
    in the Raven_errors.txt file resulting from a Raven (i.e. the C program) run.
    """

    pass


class RavenWarning(Warning):
    """
    This is a warning corresponding to a message of type "WARNING" in the Raven_errors.txt
    file resulting from a Raven (i.e. the C program) run.
    """

    pass
