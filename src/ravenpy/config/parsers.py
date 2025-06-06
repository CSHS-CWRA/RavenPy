import collections
import re
from pathlib import Path
from typing import Optional, Union

import cftime
from xarray import open_dataset

from ..config import commands as rc
from .conventions import RAVEN_OUTPUT_FMT


def parse_rv(rv: str, command: str) -> str:
    """Parse an RV file to find the value of a Command (one liner).

    Parameters
    ----------
    rv : str
        The Raven command string, e.g. "RunName".
    command : str
        Valid Raven command.

    Returns
    -------
    str:
        Command value. None if not found.
    """
    pat = re.compile(rf":{command}\s+(.+)")
    if match := re.search(pat, rv):
        return match.groups()[0]


def parse_diagnostics(fn: Path):
    """Return dictionary of performance metrics."""
    import csv

    if fn.exists():
        out = collections.defaultdict(list)

        with Path(fn).open() as f:
            reader = csv.reader(f.readlines())
            header = next(reader)
            for row in reader:
                for key, val in zip(header, row):
                    if "DIAG" in key:
                        val = float(val)
                    out[key].append(val)

            out.pop("")

        return out

    else:
        raise FileNotFoundError(fn)


def parse_nc(fn: Path):
    """Open netCDF dataset with xarray if the path is valid, otherwise return None."""
    if fn.exists():
        if fn.suffix == ".nc":
            return open_dataset(fn)
        else:
            raise NotImplementedError
    else:
        raise FileNotFoundError(fn)


def parse_solution(fn: Path, calendar: str = "PROLEPTIC_GREGORIAN"):
    """Return command objects from the model output `solution.rvc`."""
    if fn.exists():
        if fn.suffix == ".rvc":
            solution = fn.read_text()

            out = {
                "hru_state_variable_table": rc.HRUStateVariableTable.parse(solution),
                "basin_state_variables": rc.BasinStateVariables.parse(solution),
                "start_date": _time_stamp_from_solution(solution, calendar),
            }
            return out

        else:
            raise NotImplementedError
    else:
        raise FileNotFoundError(fn)


def parse_raven_messages(path):
    """Parse Raven_errors and extract the messages, structured by types."""
    messages = {
        "ERROR": [],
        "WARNING": [],
        "ADVISORY": [],
        "SIMULATION COMPLETE": False,
    }

    # The error message for an unknown command is exceptionally on two lines
    # (the second starts with a triple space)
    for m in re.findall("^([A-Z ]+) :(.+)(?:\n {3}(.+))?", path.read_text(), re.M):
        if m[0] == "SIMULATION COMPLETE":
            messages["SIMULATION COMPLETE"] = True
            continue
        msg_type = str(m[0]).strip()
        msg = f"{m[1]} {m[2]}".strip()
        if msg == "Errors found in input data. See Raven_errors.txt for details":
            # Skip this one because it's a bit circular
            continue
        messages[msg_type].append(msg)

    return messages


def output_files(run_name: str, path: Path):
    """Return path to each output file if it exists."""
    out = {}
    for k, v in RAVEN_OUTPUT_FMT.items():
        p = path / v.format(run_name=run_name + "_" if run_name is not None else "")
        if p.exists():
            out[k] = p
    return out


def parse_outputs(run_name: str, outputdir: Optional[Union[str, Path]] = None):
    """Parse outputs from model execution.

    Parameters
    ----------
    run_name : str
        RunName value identifying model outputs.
    outputdir : str or Path
        Path to model output directory. Current directory if None.

    Returns
    -------
    dict
      Dictionary holding model outputs:
        - hydrograph: xarray.Dataset
        - storage: xarray.Dataset
        - solution: Dict[str, Command]
        - diagnostics: Dict[str, list]

    Notes
    -----
    Values are set to None if no file is found.
    """
    if outputdir is None:
        outputdir = Path.cwd()
    elif isinstance(type(outputdir), str):
        outputdir = Path(outputdir)

    parser = RAVEN_OUTPUT_PARSERS

    out = {}
    files = output_files(run_name, outputdir)

    for key, path in files.items():
        if path.exists():
            if key in parser:
                out[key] = parser[key](path)
            else:
                out[key] = path.read_text()

    return out


def _time_stamp_from_solution(
    solution: str, calendar: str
) -> Optional[cftime.datetime]:
    """Return datetime from solution TimeStamp."""
    match = re.search(r":TimeStamp (\S+ \S+)", solution)
    if match:
        tt = cftime._parse_date(match.groups()[0])
        return cftime.datetime(*tt, calendar=calendar)


RAVEN_OUTPUT_PARSERS = {
    "solution": parse_solution,
    "diagnostics": parse_diagnostics,
    "storage": parse_nc,
    "hydrograph": parse_nc,
    "messages": parse_raven_messages,
}
