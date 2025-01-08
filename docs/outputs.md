---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.0
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Reading Model Outputs

ravenpy includes two relatively simples classes that expose Raven outputs as Python objects: `OutputReader` for single simulation outputs, and `EnsembleReader` to aggregate the results of multiple simulations as would be done in ensemble streamflow forecasting, for example.

(output_reader)=
## Accessing simulation results with `OutputReader`
Each Raven simulation creates output files in a directory. These outputs can be exposed as Python objects using the `ravenpy.OutputReader` class. `OutputReader` takes a `run_name` and `path` arguments, and returns an object with a number of read-only properties:
- `files`: dictionary of file path keyed by input kind;
- `solution`: dictionary with final model states for HRUs and sub-basins, as well as end date;
- `diagnostics`: dictionary storing the values of evaluation metrics;
- `hydrograph`: `xr.Dataset` of the simulated hydrograph;
- `storage`: `xr.Dataset` of the simulated storage variables;
- `messages`: the content of `Raven_errors.txt`;
- `path`: path to the output directory.

Note that `Emulator.run` returns an `OutputReader` instance, so a typical ravenpy workflow looks like this:

```{code-cell} ipython3
:tags: [skip-execution]

from ravenpy.config import emulators
from ravenpy import Emulator

# Configure model simulation
conf = emulators.HMETS(**kwds)

# Run the model
out = Emulator(conf).run()

# Look at the results
out.hydrograph.q_sim.plot()
```

Note also that the `run_name` parameter should reflect the value of the `:RunName` configuration option. If `:RunName` is not configured, then the `run_name` argument should be left to its default `None` value.

(ensemble_reader)=
## Accessing results from multiple simulations using `EnsembleReader`

Together, the `Config` and `Emulator` classes makes it fairly simple to create simulation ensembles. For example, to run the same model with different parameters, you could do something like:

```{code-cell} ipython3
:tags: [skip-execution]

# Output directory for all simulations
from pathlib import Path
p = Path("/tmp/ensemble")

# Create base model configuration
conf = HMETS(**kwds)

# Run the model for each parameter set in `params`
runs = [Emulator(conf.set_params(param), workdir=p / f"m{i}").run() for i, param in enumerate(params)]
```

Now `runs` stores a list of `OutputReader` instances. The time series stored in `hydrograph` and `storage` can be concatenated together using the `EnsembleReader` class:

```{code-cell} ipython3
:tags: [skip-execution]

from ravenpy import EnsembleReader

ens = EnsembleReader(runs=runs, dim="parameters")
ens.hydrograph.q_sim
```

where `q_sim` is a `xarray.DataArray` with dimensions `('time', 'parameters')`. An `EnsembleReader` can also be created from a list of simulation output paths:

```{code-cell} ipython3
:tags: [skip-execution]

# Create list of output paths using glob
paths = p.glob("**/output")

ens = EnsembleReader(paths=paths, dim="parameters")
```
