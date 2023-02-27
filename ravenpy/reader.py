import collections
import csv
import xarray as xr
from pathlib import Path
from config.rv import RVC


def parse(run_name, outputdir: Path):
    diag = parse_diagnostics(outputdir / f"{run_name}_Diagnostics.csv")
    sol = parse_solution(outputdir / f"{run_name}_solution.rvc")
    hydro = parse_nc(outputdir / f"{run_name}_Hydrographs.nc")
    storage = parse_nc(outputdir / f"{run_name}_WatershedStorage.nc")

    return diag, sol, hydro, storage

def parse_diagnostics(fn):
    """Return dictionary of performance metrics."""
    out = collections.defaultdict(list)

    with open(fn) as f:
        reader = csv.reader(f.readlines())
        header = next(reader)
        for row in reader:
            for (key, val) in zip(header, row):
                if "DIAG" in key:
                    val = float(val)  # type: ignore
                out[key].append(val)

        out.pop("")

    return out

def parse_nc(fn: Path):
    if fn.suffix == ".nc":
        return xr.open_dataset(fn)
    else:
        raise NotImplementedError


def parse_solution(fn: Path):
    if fn.suffix == ".rvc":
        return RVC.from_solution(fn.read_text())
    else:
        raise NotImplementedError
