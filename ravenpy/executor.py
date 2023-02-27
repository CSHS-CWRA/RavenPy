import os
import shutil
import tempfile
import subprocess
import re
from pathlib import Path
from typing import Union
from warnings import warn

RAVEN_EXEC_PATH = os.getenv("RAVENPY_RAVEN_BINARY_PATH") or shutil.which("raven")

def run(identifier: str,
        configdir: Union[str, Path],
        outputdir: Union[str, Path] = None,
        overwrite: bool = True):
    """
    Run Raven for model configuration.

    Parameters
    ----------
    identifier : str
      Configuration files name.
    configdir : str, Path
      Path to configuration files directory.
    outputdir: str, Path
      Path to model simulation output. If None, will write to configdir/output.
    """
    if not RAVEN_EXEC_PATH:
        raise RuntimeError(
            "Could not find raven binary in PATH, and RAVENPY_RAVEN_BINARY_PATH env variable is not set"
        )

    configdir = Path(configdir)
    if not configdir.exists():
        raise IOError("Workdir should include configuration files.")

    outputdir = Path(outputdir or "output")
    if not outputdir.is_absolute():
        outputdir = configdir / outputdir

    if overwrite and outputdir.exists():
        shutil.rmtree(str(outputdir))

    if not outputdir.exists():
        os.makedirs(str(outputdir))

    cmd = [RAVEN_EXEC_PATH, identifier, "-o", str(outputdir)]

    process = subprocess.Popen(
            cmd,
            cwd=configdir,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            universal_newlines=True,
        )

    process.communicate(input="\n")
    process.wait()

    messages = extract_raven_messages(configdir)

    if messages["ERROR"]:
        raise RavenError("\n".join(messages["ERROR"]))

    for msg in messages["WARNING"]:
        warn(msg, category=RavenWarning)

    assert messages["SIMULATION COMPLETE"]


def extract_raven_messages(path):
    """
    Parse all the Raven_errors and extract the messages, structured by types.
    """
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
            if (
                msg
                == "Errors found in input data. See Raven_errors.txt for details"
            ):
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

