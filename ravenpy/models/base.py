"""
Base classes
------------

The `Raven` class is the base class that implements model setup, execution and output retrieval, while the `Ostrich`
class is the base class adapting `Raven` to work with the Ostrich calibration tool.

"""
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

import numpy as np
import xarray as xr
from numpy.distutils.misc_util import is_sequence

import ravenpy
from ravenpy.config.commands import (
    DataCommand,
    GriddedForcingCommand,
    HRUsCommand,
    ObservationDataCommand,
    StationForcingCommand,
)
from ravenpy.config.rvs import RVC, Config

RAVEN_EXEC_PATH = os.getenv("RAVENPY_RAVEN_BINARY_PATH") or shutil.which("raven")
OSTRICH_EXEC_PATH = os.getenv("RAVENPY_OSTRICH_BINARY_PATH") or shutil.which("ostrich")

RAVEN_NO_DATA_VALUE = -1.2345


class RavenError(Exception):
    """
    This is an error that is meant to be raised whenever a message of type "ERROR" is found
    in the Raven_errors.txt file resulting from a Raven (i.e. the C program) run.
    """

    pass


class Raven:
    """RAVEN hydrological model wrapper.

    This class is used to run the RAVEN model from user-provided configuration files. It can also be subclassed with
    configuration templates for emulated models, allowing direct calls to the models.

    r = Raven('/tmp/testdir')
    r.configure()
    """

    _parallel_parameters = [
        "params",
        "hru_state",
        "basin_state",
        "nc_index",
        "area",
        "elevation",
        "latitude",
        "longitude",
        "region_id",
    ]

    # This is just to satisfy mypy, which wants to know about this internal class defined
    # by the emulators (which are Raven subclasses)
    Params: Any

    def __init__(
        self,
        workdir: Union[str, Path] = None,
        identifier: str = "raven-generic",
        description: str = None,
    ):
        """Initialize the RAVEN model.

        Directory for the model configuration and outputs. If None, a temporary directory will be created.
        """

        if not RAVEN_EXEC_PATH:
            raise RuntimeError(
                "Could not find raven binary in PATH, and RAVENPY_RAVEN_BINARY_PATH env variable is not set"
            )

        if not OSTRICH_EXEC_PATH:
            raise RuntimeError(
                "Could not find ostrich binary in PATH, and RAVENPY_OSTRICH_BINARY_PATH env variable is not set"
            )

        self.raven_exec = RAVEN_EXEC_PATH
        self.ostrich_exec = OSTRICH_EXEC_PATH

        # Get version from Raven binary CLI output
        out = subprocess.check_output([self.raven_exec], input="\n", text=True)
        match = re.search(r"Version (\S+) ", out)
        if match:
            self.raven_version = match.groups()[0]
        else:
            raise AttributeError(f"Raven version not found: {out}")

        self.workdir = Path(workdir or tempfile.mkdtemp())

        # Individual files for all simulations
        self.ind_outputs: Dict[str, List[Path]] = {}
        # Aggregated files
        self.outputs: Dict[str, Union[Path, str]] = {}

        # Explicit paths of every rendered RV file
        self._rv_paths: List[Path] = []

        # Directory logic
        # Top directory inside workdir. This is where Ostrich and its config and templates are stored.
        self.model_dir = "model"  # Path to the model configuration files.
        self.final_dir = "final"
        self.output_dir = "output"

        # This value is used as the stem of the generated RV file names.
        # If it's not specified by the constructor caller, it will be either:
        # (1) "raven-generic" or "ostrich-generic" for `models.base` classes
        # (2) the name in lowercase for `models.emulators` classes
        # (3) the common stem of the filenames passed to `models.base.Raven.configure`, if used
        self.identifier = identifier
        self.description = description

        self.exec_path = self.workdir / "exec"
        self.final_path = self.workdir / self.final_dir

        # Parallel simulations (one Raven config folder will be created for each)
        self._psim = 0
        self._pdim = ""  # Parallel dimension (either initparam, params or region)

        self.config = Config(model=self)

    @property
    def output_path(self):
        return self.model_path / self.output_dir

    @property
    def model_path(self):
        return self.exec_path / self.model_dir / "p{:02}".format(self.psim)

    @property
    def raven_cmd(self):
        """Path to the Raven executable."""
        return self.model_path / "raven"

    @property
    def psim(self):
        return self._psim

    @psim.setter
    def psim(self, value):
        if not isinstance(value, int):
            raise ValueError
        self.config.rvi.run_index = value
        self._psim = value

    @property
    def cmd(self):
        """This is the main executable."""
        return self.raven_cmd

    @property
    def bash_cmd(self):
        """Bash command arguments."""
        return [self.cmd, self.identifier, "-o", str(self.output_path)]

    @property
    def cmd_path(self):
        """This is the main executable."""
        return self.model_path

    def derived_parameters(self):
        """Subclassed by emulators. Defines model parameters that are a function of other parameters."""
        return

    def configure(self, fns):
        """Set configuration from existing RV files. The `self.identifier` attribute will be updated
        as the stem of the input files (which must be common to the set).
        """
        if not is_sequence(fns):
            fns = [fns]
        stems = set()
        for fn in map(Path, fns):
            identifier = fn
            while identifier.suffixes:
                identifier = Path(identifier.stem)
            stems.add(identifier)
        assert len(stems) == 1
        # Update identifier with stem
        self.identifier = str(next(iter(stems)))
        for fn in map(Path, fns):
            self.config.set_rv_file(fn)

    def _dump_rv(self):
        """Write configuration files to disk."""

        # identifier = self.config.identifier

        for rvx in ["rvt", "rvh", "rvp", "rvc", "rvi"]:
            rvo = getattr(self.config, rvx)
            if rvo.is_ostrich_tmpl:
                fn = self.exec_path / f"{self.identifier}.{rvx}.tpl"
            else:
                fn = self.model_path / f"{self.identifier}.{rvx}"
            with open(fn, "w") as f:
                self._rv_paths.append(fn)
                content = rvo.content or rvo.to_rv()
                assert (
                    content.strip()
                ), f"{rvx} has no content! (did you forget to use `RV.set_tmpl`?)"
                f.write(content)

    def setup(self, overwrite=False):
        """Create directory structure to store model input files, executable and output results.

        Model configuration files and time series inputs are stored directly in the working directory.

        workdir/  # Created by PyWPS. Is considered the model path.
        model/
        output/

        """
        if overwrite:
            if self.model_path.exists():
                shutil.rmtree(str(self.exec_path))
            if self.final_path.exists():
                shutil.rmtree(str(self.final_path))

        # Create general subdirectories
        if not self.exec_path.exists():
            os.makedirs(str(self.exec_path))  # workdir/exec
        if not self.final_path.exists():
            os.makedirs(str(self.final_path))  # workdir/final

    def setup_model_run(self, ts):
        """Create directory structure to store model input files, executable and output results.

        Parameters
        ----------
        ts : sequence
          Paths to input forcing files.
        index : int
          Run index (starts at 1)
        """
        # Compute derived parameters
        self.derived_parameters()

        # Write configuration files in model directory
        if not self.model_path.exists():
            os.makedirs(self.model_path)
            os.makedirs(self.output_path)

        self._dump_rv()

        # Create symbolic link to input files
        for fn in ts:
            if not (self.model_path / Path(fn).name).exists():
                os.symlink(str(fn), str(self.model_path / Path(fn).name))

        # Create symbolic link to Raven executable
        if not self.raven_cmd.exists():
            os.symlink(self.raven_exec, str(self.raven_cmd))

        return self.bash_cmd

    def run(self, ts, overwrite=False, **kwds):
        """Run the model.

        Parameters
        ----------
        ts : path or sequence
          Sequence of input file paths. Symbolic links to those files will be created in the model directory.
        overwrite : bool
          Whether or not to overwrite existing model and output files.
        **kwds : dict
          Raven parameters used to fill configuration file templates.

        Create a work directory with a model/ and output/ subdirectories, write the configuration files in model/ and
        launch the Raven executable. If the configuration files are templates, values can be formatted by passing
        dictionaries keyed by their extension.

        Examples
        --------
        >>> r = Raven()
        >>> r.configure(rvi='path to template', rvp='...'}
        >>> r.run(ts, start_date=dt.datetime(2000, 1, 1), area=1000, X1=67)

        """
        if isinstance(ts, (str, Path)):
            ts = [ts]

        # Support legacy interface for single HRU emulator
        hru_attrs = {}
        for k in ["area", "latitude", "longitude", "elevation"]:
            v = kwds.pop(k, None)
            if v:
                # It seems that `v` is a list when running via a WPS interface
                hru_attrs[k] = v[0] if isinstance(v, list) else v
        if hru_attrs:
            assert len(self.config.rvh.hrus) == 1
            self.config.rvh.hrus = (replace(self.config.rvh.hrus[0], **hru_attrs),)

        # Case for potentially parallel parameters
        pdict = {}
        for p in self._parallel_parameters:
            val = kwds.pop(p, None)
            if val is not None and p == "params":
                assert hasattr(self, "Params")  # make sure we are in an emulator
                # Special case where we have `Params(..)` or `[Params(), ..]`
                lval = [val] if not is_sequence(val) else val
                if isinstance(lval[0], self.Params):
                    pdict[p] = np.atleast_1d(val)
                else:
                    pdict[p] = np.atleast_2d(val)
            else:
                pdict[p] = np.atleast_1d(val)

        # Number of parallel loops is dictated by the number of parallel parameters or nc_index.
        plen = {pp: len(pdict[pp]) for pp in self._parallel_parameters + ["nc_index"]}

        # Find the longest parallel array and its length
        longer, nloops = max(plen.items(), key=operator.itemgetter(1))

        # Assign the name of the parallel dimension
        # nbasins is set by RavenC++
        if nloops > 1:
            self._pdim = {
                "params": "params",
                "hru_state": "state",
                "basin_state": "state",
                "nc_index": "nbasins",
            }[longer]

        for key, val in pdict.items():
            if len(val) not in [1, nloops]:
                raise ValueError(
                    "Parameter {} has incompatible dimension: {}. "
                    "Should be 1 or {}.".format(key, len(val), nloops)
                )

        # Resize parallel parameters to the largest size
        for key, val in pdict.items():
            if len(val) == 1:
                pdict[key] = val.repeat(nloops, axis=0)

        # Use rvc file to set model state, if any
        rvc = kwds.pop("rvc", None)
        if rvc:
            self.resume(solution=rvc)

        # Update non-parallel parameter objects
        for key, val in kwds.items():
            self.config.update(key, val)

        ts_ncs = [f for f in ts if Path(f).suffix.startswith(".nc")]

        if ts_ncs:
            self.config.rvi.configure_from_nc_data(ts_ncs)

        if ts_ncs and self.config.rvt._auto_nc_configure:
            self.config.rvt.configure_from_nc_data(ts_ncs)

        # Loop over parallel parameters - sets self.rvi.run_index
        procs = []
        for self.psim in range(nloops):
            for key, val in pdict.items():
                if val[self.psim] is not None:
                    if key == "hru_state":
                        self.config.rvc.set_hru_state(val[self.psim])
                    elif key == "basin_state":
                        self.config.rvc.set_basin_state(val[self.psim])
                    else:
                        self.config.update(key, val[self.psim])

            cmd = self.setup_model_run(tuple(map(Path, ts)))

            procs.append(
                subprocess.Popen(
                    cmd,
                    cwd=self.cmd_path,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    universal_newlines=True,
                )
            )

        return procs

    def __call__(self, ts, overwrite=False, **kwds):
        self.setup(overwrite)

        procs = self.run(ts, overwrite, **kwds)

        for proc in procs:
            # When Raven errors right away (for instance if it's missing an RV file)
            # it asks for a RETURN to exit
            proc.communicate(input="\n")
            proc.wait()

        messages = self.extract_raven_messages()

        if messages["ERROR"]:
            raise RavenError("\n".join(messages["ERROR"]))

        assert messages["SIMULATION COMPLETE"]

        self.parse_results()

    def resume(self, solution=None):
        """Set the initial state to the state at the end of the last run.

        Parameters
        ----------
        solution : str, Path
          Path to solution file. If None, will use solution from last model run if any.
        """
        if solution is None:
            fn = self.outputs["solution"]
        else:
            fn = solution

        self.config.rvc.parse_solution(Path(fn).read_text())

    def parse_results(self, path=None, run_name=None):
        """Store output files in the self.outputs dictionary."""
        # Output files default names. The actual output file names will be composed of the run_name and the default
        # name.
        path = path or self.exec_path
        run_name = run_name or self.config.rvi.run_name or ""
        patterns = {
            "hydrograph": f"{run_name}*Hydrographs.nc",
            "storage": f"{run_name}*WatershedStorage.nc",
            "solution": f"{run_name}*solution.rvc",
            "diagnostics": f"{run_name}*Diagnostics.csv",
        }

        for key, pattern in patterns.items():
            # There are no diagnostics if a streamflow time series is not provided.
            try:
                fns = self._get_output(pattern, path=path)
            except UserWarning as exc:
                if key != "diagnostics":
                    raise exc
                else:
                    continue

            fns.sort()
            self.ind_outputs[key] = fns
            self.outputs[key] = self._merge_output(fns, pattern[1:])

        self.outputs["rv_config"] = self._merge_output(self._rv_paths, "rv.zip")

    def _merge_output(self, files, name):
        """Merge multiple output files into one if possible, otherwise return a list of files."""
        # If there is only one file, return its name directly.
        from .multimodel import RavenMultiModel

        if len(files) == 1:
            return files[0]

        # Otherwise try to create a new file aggregating all files.
        outfn = self.final_path / name

        if name.endswith(".nc") and not isinstance(self, RavenMultiModel):
            ds = [xr.open_dataset(fn) for fn in files]
            try:
                # We aggregate along the pdim dimensions.
                out = xr.concat(ds, self._pdim, data_vars="all")
                out.to_netcdf(outfn)
                return outfn
            except (ValueError, KeyError):
                pass

        # Let's zip the files that could not be merged.
        outfn = outfn.with_suffix(".zip")

        # Find the lower file parts level at which there are differences among files.
        i = get_diff_level(files)

        # Try to create a zip file
        with zipfile.ZipFile(outfn, "w") as f:
            for fn in files:
                len(fn.parts)
                f.write(fn, arcname=fn.relative_to(Path(*fn.parts[:i])))

        return outfn

    def extract_raven_messages(self):
        """
        Parse all the Raven_errors and extract the messages, structured by types.
        """
        err_filepaths = self.exec_path.rglob("Raven_errors.txt")
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
                if (
                    msg
                    == "Errors found in input data. See Raven_errors.txt for details"
                ):
                    # Skip this one because it's a bit circular
                    continue
                messages[msg_type].append(msg)  # type: ignore

        return messages

    def _get_output(self, pattern, path):
        """Match actual output files to known expected files.

        Return a dictionary of file paths for each expected input.
        """
        files = list(path.rglob(pattern))

        if len(files) == 0:
            if not self.config.rvi.suppress_output:
                raise UserWarning("No output files for {} in {}.".format(pattern, path))

        return [f.absolute() for f in files]

    @property
    def q_sim(self):
        """Return a view of the hydrograph time series.

        This view will be overwritten by successive calls to `run`. To make a copy of this DataArray that will
        persist in memory, use `q_sim.copy(deep=True)`.
        """
        if isinstance(self.hydrograph, list):
            return [h.q_sim for h in self.hydrograph]

        return self.hydrograph.q_sim

    @property
    def hydrograph(self):
        """Return a view of the current output file.

        If the model is run multiple times, hydrograph will point to the latest version. To store the results of
        multiple runs, either create different model instances or explicitly copy the file to another disk location.
        """
        hydrograph = cast(Path, self.outputs["hydrograph"])
        if hydrograph.suffix == ".nc":
            return xr.open_dataset(hydrograph)
        elif hydrograph.suffix == ".zip":
            return [xr.open_dataset(fn) for fn in self.ind_outputs["hydrograph"]]
        else:
            raise ValueError

    @property
    def storage(self):
        storage = cast(Path, self.outputs["storage"])
        if storage.suffix == ".nc":
            return xr.open_dataset(storage)
        elif storage.suffix == ".zip":
            return [xr.open_dataset(fn) for fn in self.ind_outputs["storage"]]
        else:
            raise ValueError

    @property
    def solution(self):
        solution = cast(Path, self.outputs["solution"])
        if solution.suffix == ".rvc":
            return RVC.create_solution(solution.read_text())
        elif solution.suffix == ".zip":
            return [
                RVC.create_solution(fn.read_text())
                for fn in self.ind_outputs["solution"]
            ]

    def get_final_state(self, hru_index=1, basin_index=1):
        """Return model state at the end of simulation.

        Parameters
        ----------
        hru_index : None, int
          Set index value or None to get all HRUs.
        basin_index : None, int
          Set index value or None to get all basin states.
        """
        solution = self.solution
        if isinstance(solution, RVC):
            return solution.hru_states[hru_index], solution.basin_states[basin_index]
        else:
            states = [
                (sol.hru_states[hru_index], sol.basin_states[basin_index])
                for sol in solution
            ]
            return zip(*states)

    @property
    def diagnostics(self):
        """Return a nested dictionary of performance metrics keyed by diagnostic name and period. The default period
        is called "ALL".
        """
        diag = []
        out = collections.defaultdict(list)
        for fn in self.ind_outputs["diagnostics"]:
            with open(fn) as f:
                reader = csv.reader(f.readlines())
                header = next(reader)
                for row in reader:
                    for (key, val) in zip(header, row):
                        if "DIAG" in key:
                            val = float(val)  # type: ignore
                        out[key].append(val)

                out.pop("")

            diag.append(out)
        return diag if len(diag) > 1 else diag[0]


class Ostrich(Raven):
    """Wrapper for OSTRICH calibration of RAVEN hydrological model.

    This class is used to calibrate RAVEN model using OSTRICH from user-provided configuration files. It can also be
    subclassed with configuration templates for emulated models, allowing direct calls to the models.

    Parameters
    ----------
    conf:
      The rv configuration files + Ostrict ostIn.txt.
    tpl:
      The Ostrich templates.

    Examples
    --------
    >>> r = Ostrich('/tmp/testdir')
    >>> r.configure()
    """

    def __init__(self, *args, **kwds):
        kwds["identifier"] = kwds.get("identifier", "ostrich-generic")
        super().__init__(*args, **kwds)

    @property
    def model_path(self):
        return self.exec_path / self.model_dir

    @property
    def ostrich_cmd(self):
        """OSTRICH executable path."""
        return self.exec_path / "ostrich"

    @property
    def cmd(self):
        """OSTRICH executable path."""
        return self.ostrich_cmd

    @property
    def cmd_path(self):
        """This is the main executable."""
        return self.exec_path

    @property
    def proc_path(self):
        """Path to Ostrich parallel process directory."""
        return self.exec_path / "processor_0"  # /'model' / 'output' ?

    def write_save_best(self):
        fn = self.exec_path / "save_best.sh"
        fn.write_text(save_best)
        make_executable(fn)

    def write_ostrich_runs_raven(self):
        fn = self.exec_path / "ostrich-runs-raven.sh"
        fn.write_text(ostrich_runs_raven.format(identifier=self.identifier))
        make_executable(fn)

    def setup(self, overwrite=False):
        """Create directory structure to store model input files, executable, and output results.

        Model configuration files and time series inputs are stored directly in the working directory.
        At each Ostrich loop, configuration files (original and created from templates are copied into model).
        """
        Raven.setup(self, overwrite)

        os.makedirs(str(self.final_path), exist_ok=True)

        self.write_ostrich_runs_raven()
        self.write_save_best()

        # Create symbolic link to executable
        if not self.cmd.exists():
            os.symlink(self.ostrich_exec, str(self.cmd))

    def configure(self, fns):
        """Set configuration from existing RV files. The `self.identifier` attribute will be updated
        as the stem of the input files (which must be common to the set).
        """
        if not is_sequence(fns):
            fns = [fns]
        stems = set()
        for fn in map(Path, fns):
            identifier = fn
            while identifier.suffixes:
                identifier = Path(identifier.stem)
            # OST-related files have a different stem
            if not str(identifier).lower().startswith("ost"):
                stems.add(identifier)
        if not stems:
            # Special case
            assert fns[0].name == "OstRandomNumbers.txt"
        else:
            assert len(stems) == 1
            # Update identifier with stem
            self.identifier = str(next(iter(stems)))
        for fn in map(Path, fns):
            self.config.set_rv_file(fn)

    def _dump_rv(self):
        """write configuration files to disk."""

        super()._dump_rv()

        # ostIn.txt
        fn = self.exec_path / "ostIn.txt"
        with open(fn, "w") as f:
            self._rv_paths.append(fn)
            content = self.config.ost.content or self.config.ost.to_rv()
            assert (
                content.strip()
            ), "OST has no content! (did you forget to use `RV.set_tmpl`?)"
            f.write(content)

        # OstRandomNumbers.txt
        if self.config.ost.random_numbers_path:
            fn = self.exec_path / "OstRandomNumbers.txt"
            with open(fn, "w") as f:
                f.write(self.config.ost.random_numbers_path.read_text())
            self._rv_paths.append(fn)

    def parse_results(self):
        """Store output files in the self.outputs dictionary."""
        # Output files default names. The actual output file names will be composed of the run_name and the default
        # name.
        Raven.parse_results(self, path=self.final_path)

        patterns = {
            "params_seq": "OstModel?.txt",
            "calibration": "OstOutput?.txt",
        }

        # Store output file names in dict
        for key, pattern in patterns.items():
            fns = self._get_output(pattern, path=self.exec_path)
            if len(fns) == 1:
                fns = fns[0]
            self.outputs[key] = fns

        try:
            if is_dataclass(self.calibrated_params):
                self.outputs["calibparams"] = ", ".join(
                    map(str, astuple(self.calibrated_params))
                )
            else:
                self.outputs["calibparams"] = ", ".join(
                    map(str, self.calibrated_params)
                )
        except (AttributeError, TypeError):
            err = self.parse_errors()
            raise UserWarning(err)

    def parse_errors(self):
        try:
            raven_err = self._get_output("OstExeOut.txt", path=self.exec_path)[
                0
            ].read_text()
        except UserWarning:  # Read in processor_0 directory instead.
            try:
                raven_err = self._get_output("OstExeOut.txt", path=self.proc_path)[
                    0
                ].read_text()
            except UserWarning:
                raven_err = ""

        try:
            ost_err = self._get_output("OstErrors?.txt", path=self.exec_path)[
                0
            ].read_text()
        except UserWarning:  # Read in processor_0 directory instead.
            ost_err = self._get_output("OstErrors?.txt", path=self.proc_path)[
                0
            ].read_text()

        return f"{ost_err}\n{raven_err}"

    def parse_optimal_parameter_set(self):
        """Return dictionary of optimal parameter set."""
        txt = open(self.outputs["calibration"]).read()
        ops = re.search(r".*Optimal Parameter Set(.*?)\n{2}", txt, re.DOTALL).groups()  # type: ignore
        p = re.findall(r"(\w+)\s*:\s*([\S]+)", ops[0])
        return OrderedDict((k, float(v)) for k, v in p)

    def ost2raven(self, ostrich_params):
        """Return model parameters.

        Parameters
        ----------
        ostrich_params: dict

        If this method is used in the context of an emulator subclass, it will return
        an instance of the Params dataclass of the emulator (it will also make sure that the param
        name conversion from Ostrich to Raven is performed); if not, it will return the list of values.

        """

        if hasattr(self, "Params"):
            # We are in an emulator subclass, so a properly typed Params dataclass is available
            raven_params = {}
            # Perform Ostrich to Raven param name conversion if needed
            o2r = getattr(self, "ostrich_to_raven_param_conversion", {})
            r2o = {r: o for o, r in o2r.items()}
            for f in fields(self.Params):
                raven_params[f.name] = ostrich_params[r2o.get(f.name, f.name)]
            return self.Params(**raven_params)
        else:
            # We are using generic Ostrich
            return ostrich_params.values()

    @property
    def calibrated_params(self):
        """The dictionary of optimal parameters estimated by Ostrich."""
        ops = self.parse_optimal_parameter_set()
        return self.ost2raven(ops)

    @property
    def obj_func(self):
        return np.loadtxt(self.outputs["params_seq"], skiprows=1)[-1, 1]

    @property
    def optimized_parameters(self):
        """These are the raw parameters returned by Ostrich."""
        return np.loadtxt(self.outputs["params_seq"], skiprows=1)[-1, 2:]


def get_diff_level(files):
    """Return the lowest hierarchical file parts level at which there are differences among file paths."""

    for i, parts in enumerate(zip(*[f.parts for f in files])):
        if len(set(parts)) > 1:
            return i


def make_executable(fn):
    """Make file executable."""
    st = os.stat(fn)
    os.chmod(fn, st.st_mode | stat.S_IEXEC)


def get_average_annual_runoff(
    nc_file_path,
    area_in_m2,
    time_dim="time",
    obs_var="qobs",
    na_value=RAVEN_NO_DATA_VALUE,
):
    """
    Compute the average annual runoff from observed data.
    """
    with xr.open_dataset(nc_file_path) as ds:
        qobs = ds.where(ds[obs_var] != na_value)[obs_var]
        qobs *= 86400.0  # convert m**3/s to m**3/d
        axis = qobs.dims.index(time_dim)
        # avg daily runoff [m3/d] for each year in record
        qyear = np.nanmean(qobs.groupby("time.year").mean("time"), axis=axis)
        qyear = qyear / area_in_m2 * 365 * 1000.0  # [mm/yr] for each year in record
        qyear = np.mean(qyear)  # [mm/yr] mean over all years in record

    return qyear


# TODO: Configure this according to the model_path and output_path.
save_best = """#!/bin/bash

set -e

cp ./model/*.rv?  ../../final/
cp ./model/output/* ../../final/

exit 0
"""

# TODO: Configure this according to raven_cmd, name and output_path.
ostrich_runs_raven = """
#!/bin/bash

set -e

cp ./*.rv? model/

./model/raven ./model/{identifier} -o ./model/output/

exit 0
"""
