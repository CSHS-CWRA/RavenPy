"""Main module."""

import collections
import os
import shutil
import subprocess
from pathlib import Path
from typing import Union
from warnings import warn

import xarray as xr

from ravenpy.new_config import conventions, parsers
from ravenpy.new_config.rvs import Config

RAVEN_EXEC_PATH = os.getenv("RAVENPY_RAVEN_BINARY_PATH") or shutil.which("raven")


class Emulator:
    def __init__(
        self, config: Config, workdir: Union[Path, str], modelname=None, overwrite=False
    ):
        """
        Convenience class to work with the Raven modeling framework.

        Parameters
        ----------
        config : Config
          Emulator Config instance fully parameterized, i.e. without symbolic expressions.
        path : str, Path
          Path to rv files and model outputs.
        modelname : str, None
          File name stem of configuration files: `<modelname>.rv*`.
        """
        self._config = config.copy()
        self._workdir = Path(workdir)
        self._output = None  # Model output path
        self._modelname = modelname
        self._rv = None  # Model config files
        self._overwrite = overwrite

        # Emulator config is made immutable to avoid complex behavior.
        self._config.__config__.allow_mutation = False

    def write_rv(self, overwrite=False):
        """Write the configuration files to disk."""
        if self._overwrite is not None:
            overwrite = self._overwrite

        self._rv = self._config.write_rv(
            workdir=self.workdir, modelname=self.modelname, overwrite=overwrite
        )
        if self.modelname is None:
            # modelname default set by `write_rv`
            self._modelname = self._rv["rvi"].stem
        return self._rv

    def run(self) -> "OutputReader":
        """Run the model. This will write RV files if not already done."""
        if not (self.workdir / f"{self.modelname}.rvi").exists():
            self.write_rv()

        self._output_path = run(self.modelname, self.workdir, "output")
        self._output = OutputReader(self.config.run_name, path=self._output_path)
        return self._output

    @property
    def config(self) -> Config:
        """Read-only model configuration."""
        return self._config

    @property
    def workdir(self) -> Path:
        """Path to RV files and output sub-directory."""
        return self._workdir

    @property
    def output_path(self) -> Path:
        """Path to model outputs."""
        return self._output_path

    @property
    def modelname(self) -> str:
        """File name stem of configuration files."""
        return self._modelname

    @property
    def output(self) -> "OutputReader":
        """Return simulation output object."""
        return self._output

    def resume(self) -> Config:
        """Return new model configuration using state variables from the end of the run."""
        out = self.config.set_solution(self._output.files["solution"])
        out.__config__.allow_mutation = True
        return out


class OutputReader:
    def __init__(self, run_name: str = None, path: Path = None):
        """Class facilitating access to Raven model output.

        Parameters
        ----------
        run_name : str, None
          Simulation name, if any is specified by the `RunName` configuration.
        path : str, Path
          Output directory where model results are stored. Defaults to the current directory.
        """
        self._run_name = run_name
        self._path = Path(path) if path else Path.cwd()

    @property
    def files(self) -> dict:
        """Return paths to output files."""
        return parsers.output_files(self._run_name, self._path)

    @property
    def objs(self) -> dict:
        """Return model outputs as Python objects.

        - diagnostics: dict
        - hydrograph: xr.Dataset
        - storage: xr.Dataset
        - solution:

        """
        return parsers.parse_outputs(self._run_name, self._path)

    @property
    def solution(self) -> str:
        """Return solution file content."""
        solution = self.files.get("solution")
        if solution:
            return parsers.parse_solution(solution)

    @property
    def diagnostics(self) -> dict:
        """Return model diagnostics."""
        diag = self.files.get("diagnostics")
        if diag:
            return parsers.parse_diagnostics(diag)

    @property
    def hydrograph(self) -> xr.Dataset:
        """Return the hydrograph."""
        h = self.files.get("hydrograph")
        if h:
            return parsers.parse_nc(h)

    @property
    def storage(self) -> xr.Dataset:
        s = self.files.get("storage")
        if s:
            return parsers.parse_nc(s)

    @property
    def messages(self) -> str:
        msg = self.files.get("messages")
        if msg:
            return msg.read_text()

    @property
    def path(self) -> Path:
        """Path to output directory."""
        return self._path


class EnsembleReader:
    def __init__(self, run_name, paths, dim="member"):
        """
        Class facilitating access to ensemble of Raven outputs.

        Parameters
        ----------
        run_name: str
          Name of simulation.
        paths: list
          List of output paths.
        """
        self._paths = paths
        self._outputs = [OutputReader(run_name, p) for p in paths]
        self._dim = dim

    @property
    def files(self):
        out = collections.defaultdict(list)
        for o in self._outputs:
            for k, v in o.files.items():
                out[k].append(v)
        return out

    @property
    def storage(self):
        return xr.concat(
            [xr.open_dataset(f) for f in self.files["storage"]], dim=self._dim
        )

    @property
    def hydrograph(self):
        return xr.concat(
            [xr.open_dataset(f) for f in self.files["hydrograph"]], dim=self._dim
        )


def run(
    modelname: str,
    configdir: Union[str, Path],
    outputdir: Union[str, Path] = None,
    overwrite: bool = True,
):
    """
    Run Raven given the path to an existing model configuration.

    Parameters
    ----------
    modelname : str
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
    cmd = [RAVEN_EXEC_PATH, modelname, "-o", str(outputdir)]

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
    messages = parsers.parse_raven_messages(outputdir / "Raven_errors.txt")

    for msg in messages["WARNING"] + messages["ADVISORY"]:
        warn(msg, category=RavenWarning)

    if messages["ERROR"]:
        raise RavenError("\n".join(messages["ERROR"]))

    return outputdir


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
