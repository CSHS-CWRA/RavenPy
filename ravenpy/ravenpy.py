"""Main module."""

import collections
import os
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path
from typing import List, Optional, Union
from warnings import warn

import xarray as xr

from .config import parsers
from .config.rvs import Config

RAVEN_EXEC_PATH = os.getenv("RAVENPY_RAVEN_BINARY_PATH") or shutil.which("raven")


class Emulator:
    def __init__(
        self,
        config: Config,
        workdir: Optional[Union[str, os.PathLike]] = None,
        modelname: Optional[str] = None,
        overwrite: bool = False,
    ):
        """Convenience class to work with the Raven modeling framework.

        Parameters
        ----------
        config : Config
            Emulator Config instance fully parameterized, i.e. without symbolic expressions.
        workdir : Path or str
            Path to rv files and model outputs. If None, create a temporary directory.
        modelname : str, None
            File name stem of configuration files: `<modelname>.rv*`.
        overwrite : bool
            If True, overwrite existing files.
        """
        self._config = config.copy(deep=True)
        self._workdir = Path(workdir or tempfile.mkdtemp())
        self._modelname = modelname
        self.overwrite = overwrite
        self._output = None  # Model output path
        self._output_path = None

        # Write model config files
        self._rv = self._config.write_rv(
            workdir=self.workdir, modelname=self.modelname, overwrite=overwrite
        )
        # Grab modelname in case it was set by config.write_rv
        self._modelname = self._rv["rvi"].stem

    def run(self, overwrite: bool = False) -> "OutputReader":
        """Run the model. This will write RV files if not already done.

        Parameters
        ----------
        overwrite : bool
            If True, overwrite existing files.
        """
        if not (self.workdir / f"{self.modelname}.rvi").exists():
            # FIXME: No attribute 'write_rv' on Emulator [attribute-error]
            self.write_rv(overwrite=overwrite)

        self._output_path = run(
            self.modelname, self.workdir, "output", overwrite=overwrite
        )
        self._output = OutputReader(self.config.run_name, path=self._output_path)
        return self._output

    @property
    def config(self) -> Config:
        """Read-only model configuration."""
        return self._config

    @property
    def workdir(self) -> Path:
        """Path to RV files and output subdirectory."""
        return self._workdir

    @property
    def output_path(self) -> Optional[Path]:
        """Path to model outputs."""
        if self._output_path is not None:
            return self._output_path
        warnings.warn("`output_path` not set. Model must be run first")

    @property
    def modelname(self) -> str:
        """File name stem of configuration files."""
        return self._modelname

    @property
    def output(self) -> "OutputReader":
        """Return simulation output object."""
        return self._output

    def resume(self, timestamp: bool = True) -> Config:
        """Return new model configuration using state variables from the end of the run.

        timestamp: bool
          If False, ignore time stamp information in the solution. If True, the solution
          will set StartDate to the solution's timestamp.
        """
        return self.config.set_solution(
            self.output.files["solution"], timestamp=timestamp
        )


class OutputReader:
    def __init__(
        self, run_name: Optional[str] = None, path: Optional[Union[str, Path]] = None
    ):
        """Class facilitating access to Raven model output.

        Parameters
        ----------
        run_name : str, optional
            Simulation name, if any is specified by the `RunName` configuration.
        path : str or Path
            Output directory where model results are stored. Defaults to the current directory.
        """
        self._run_name = run_name
        self._path = Path(path) if path else Path.cwd()
        self._files = parsers.output_files(self._run_name, self._path)

        # TODO: Check if no files are found. Otherwise we get cryptic errors.
        # if self._files["hydrograph"]

    @property
    def files(self) -> dict:
        """Return paths to output files."""
        return self._files

    @property
    def solution(self) -> Optional[dict]:
        """Return solution file content."""
        solution = self.files.get("solution")
        if solution:
            return parsers.parse_solution(solution)

    @property
    def diagnostics(self) -> Optional[dict]:
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
    def messages(self) -> Optional[str]:
        msg = self.files.get("messages")
        if msg:
            return msg.read_text()

    @property
    def path(self) -> Path:
        """Path to output directory."""
        return self._path


class EnsembleReader:
    def __init__(
        self,
        *,
        run_name: Optional[str] = None,
        paths: Optional[List[Union[str, os.PathLike]]] = None,
        runs: List[OutputReader] = None,
        dim: str = "member",
    ):
        """
        Class facilitating access to ensemble of Raven outputs.

        Parameters
        ----------
        run_name : str, None
            Name given to simulation, if any.
        paths : list[str or Path], optional
            List of output paths. Defaults to all directories in current directory.
        runs : List[OutputReader]
            List of OutputReader instances.
        dim : str
            Name of concatenation dimension.
        """
        self._dim = dim
        if runs is not None:
            self._outputs = runs
            self._paths = [o.path for o in self._outputs]

        else:
            if paths is None:
                paths = [p for p in Path.cwd().iterdir() if p.is_dir()]
            self._paths = [Path(p) for p in paths]
            self._outputs = [OutputReader(run_name, p) for p in self._paths]

    @property
    def files(self):
        out = collections.defaultdict(list)
        for o in self._outputs:
            for k, v in o.files.items():
                out[k].append(v)
        return out

    @property
    def storage(self):
        if len(self.files["storage"]) == 0:
            raise ValueError(
                "No file found, make sure you have the right `run_name` and output `paths`."
            )

        return xr.concat(
            [xr.open_dataset(f) for f in self.files["storage"]], dim=self._dim
        )

    @property
    def hydrograph(self):
        if len(self.files["hydrograph"]) == 0:
            raise ValueError(
                "No file found, make sure you have the right `run_name` and output `paths`."
            )
        return xr.concat(
            [xr.open_dataset(f) for f in self.files["hydrograph"]],
            dim=self._dim,
            coords="different",
        )


def run(
    modelname: str,
    configdir: Union[str, Path],
    outputdir: Optional[Union[str, Path]] = None,
    overwrite: bool = True,
    verbose: bool = False,
) -> Path:
    """Run Raven given the path to an existing model configuration.

    Parameters
    ----------
    modelname : str
        Configuration files stem, i.e. the file name without extension.
    configdir : Path or str
        Path to configuration files directory.
    outputdir : Path or str, optional
        Path to model simulation output.
        If None, will write to configdir/output.
    overwrite : bool
        If True, overwrite existing files.
    verbose : bool
        If True, always display Raven warnings. If False, warnings will only be printed if an error occurs.

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

    # Deal with errors and warnings
    messages = parsers.parse_raven_messages(outputdir / "Raven_errors.txt")

    if messages["ERROR"] or verbose:
        for msg in messages["WARNING"] + messages["ADVISORY"]:
            warn(msg, category=RavenWarning)

    if messages["ERROR"]:
        raise RavenError(
            "\n".join([f"Config directory: {configdir}"] + messages["ERROR"])
        )

    if returncode != 0:
        raise OSError(f"Raven segfaulted : \n{stdout}")

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
