---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Model Execution

ravenpy provides two mechanisms to execute a model configuration. The first one described is the `ravenpy.run` function, which expects to find a directory with RV files and returns the path to the output directory. The second one is the `ravenpy.Emulator` class which exposes more functionalities and deserves some explanations.

(emulator)=
## The Emulator class

The `Emulator` class minimally expects a fully configured `Config` instance, which it immediately writes to disk. Other optional arguments are `workdir`, the path to the work directory where configuration files are written, and `modelname`, the string used to create configuration file names. By default, `workdir` will be set to a temporary directory, so make sure to specify another location if you want to hold on to the results.

When the `run` method of an `Emulator` instance is called, it launches the model which writes simulations results in the `output` folder, then returns an instance of the `OutputReader` class. Beyond the `write_rv` and `run` methods, Emulator has a number of read-only properties:
- `config`: the model configuration `Config` instance;
- `workdir`: the path to the work directory where configuration files are written;
- `output_path`: the path to simulation results;
- `modelname`: the name given to configuration files written to disk;
- `output`: the `OutputReader` instance created out of the simulation results.

The `Emulator` has another method called `resume`, which returns a new `Config` instance where the initial states for HRUs and sub-basins are set from the `solution.rvc` file storing the states at the end of the run. Note that by default, the `StartDate` of this new configuration will be set to the date following the end of the run. Set `timestamp=False` to ignore the time stamp information of the solution file.
