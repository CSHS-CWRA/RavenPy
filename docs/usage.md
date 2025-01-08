# Usage

RavenPy is designed to make Raven easier to use from an interactive Python programming environment. It can:
- execute Raven on an existing model configuration;
- expose simulation results from one or more simulations as `xarray.Datasets`;
- create, modify and write model configurations.

In particular, RavenPy includes eight pre-configured model *emulators*.


## Running Raven from existing configuration files

To run Raven using existing configuration files (`.rv*`), simply call the `run` function with the name of the configuration file and the path to the directory storing the RV files:

```python
from ravenpy import run

output_path = run(modelname, configdir)
```

`run` simply returns the directory storing the model outputs. If Raven emits warnings, those will be printed in the console. If Raven raises errors, `run` will halt the execution with a `RavenError` and display the error messages in `Raven_errors.txt`.


## Exposing model outputs as Python objects

The model outputs can be read with the `OutputReader` class:

```python
from ravenpy import OutputReader

out = OutputReader(run_name, path=output_path)
out.hydrograph.q_sim
out.storage["Ponded Water"]
```

Note that this works only if simulated variables are stored as netCDF files, that is, if the `rvi` file includes a `:WriteNetCDFFormat` command. For more details, see {ref}`output_reader`.

The class `EnsembleReader` does the same for an ensemble of model outputs, concatenating netCDF outputs along a new dimension:

```python
from ravenpy import EnsembleReader

out = EnsembleReader(
    run_name,
    paths=[
        output_path,
    ],
    dim="ensemble_dim",
)
out.hydrograph.q_sim
out.storage["Ponded Water"]
```
For more info, see {ref}`ensemble_reader`.

## Configuring emulators

Ravenpy comes packaged with pre-configured emulators, that is, Raven model configurations that can be modified on the fly. These emulators are made out of symbolic expressions, connecting model parameters to properties and coefficients. For example, the code below creates a model configuration for emulated model GR4JCN using the parameters given, as well as a `Gauge` configuration inferred by inspecting the `meteo.nc` file.

```python
from ravenpy.config import emulators, commands

gr4jcn = emulators.GR4JCN(
    params=[0.5, -3.0, 400, 1.0, 17, 0.9], Gauge=[commands.Gauge.from_nc("meteo.nc")]
)
```

Note that `Gauge.from_nc` will only find the required information if the netCDF file complies with the `CF-Convention`. Otherwise, additional parameters have to be provided to complete the configuration.

Ravenpy includes a suite of eight emulators:
  - GR4JCN (GR4J-Cemaneige)
  - HMETS
  - HBVEC
  - Mohyse
  - Blended
  - CanadianShield
  - SACSMA
  - HYPR


## Running an emulator

The RV files for the emulator above can be inspected using the `rvi`, `rvh`, `rvp`, `rvc` and `rvt` properties, e.g. `print(gr4jcn.rvt)` will show the `rvt` file as it would be written to disk. Configuration files can then be written to disk using `gr4jcn.write_rv(workdir, modelname)`, and the model launched using the `run` function introduced before.

For convenience, `ravenpy` also proposes the `Emulator` class, designed to streamline the execution of the model and the retrieval of the results.

```python
from ravenpy import Emulator

e = Emulator(config=gr4jcn, workdir="/tmp/gr4jcn/run_1")
out = e.run()
out.hydrograph.q_sim
```
If no `workdir` is given, a temporary directory will be created, available from  The `Emulator.workdir` property. `Emulator` also has `resume` method that returns a copy of the original configuration with the internal states and start date set to the values stored in the `solution.rvc` file, which can then be used to launch another simulation following the first one. For more on this, see {ref}`emulator`.

For more information on model configuration, see {ref}`conf`.
