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

# Model configuration


Raven lets hydrologist customize how models are defined, and provides examples of model configuration that can reproduce outputs from known hydrological models: HBV-EC, HMETS, Sacramento-SMA, etc. RavenPy can be used to create Python objects representing those emulators, facilitating the setup of complex modeling experiments. The following describes the basics of emulator configuration.

## The `Config` class

RavenPy writes configuration files using the `Config` class. `Config` recognizes most Raven commands, and has mechanisms to parse and validate the inputs, and then write them in the appropriate RV file. For example:

```{code-cell} ipython3
from ravenpy.new_config import Config
conf = Config(StartDate="2023-03-31")
print(conf.rvi)
```

`Config` has recognized `StartDate` as a Raven command, and knows it should appear in the `rvi` file as a line starting with `:StartDate` followed by a date in ISO format. `StartDate` could equally have been given as a `datetime.date` or `datetime.datetime` object, and `Config` would have parsed it correctly.

Many other Raven commands are known to `Config`. To see what commands are supported, and the type of inputs they expect, consult the class docstring with `help(Config)`. All these commands are set to a default value of `None`, and won't appear in RV files unless they are given actual values. If you plan on using the `OutputReader`, make sure you set `WriteNetcdfFormat` to True.

Many Raven commands take values as strings, floats, integers or dates, but some accept lists of arguments, such as `:EvaluationMetrics`:

```{code-cell} ipython3
conf = Config(EvaluationMetrics=["NASH_SUTCLIFFE", "RMSE"])
print(conf.rvi)
```


Some commands require more complex structures, for example, the configuration for  `CustomOutput` is given as a dictionary with arguments `time_per`, `stat`, `variable` and `space_agg`:

```{code-cell} ipython3
conf = Config(CustomOutput={"time_per": "MONTHLY", "stat": "AVERAGE", "variable": "SNOW", "space_agg": "BY_BASIN"})
print(conf.rvi)
```

All Raven Commands with inputs more complex than single values or simple lists   are defined in `ravenpy.config.commands`. Consult the docstring to find out how each should be instantiated. Attributes can be given as dictionaries that `Config` will parse, as above, or as `Command` instances:

```{code-cell} ipython3
from ravenpy.new_config import commands as rc

co = rc.CustomOutput(time_per="MONTHLY", stat="AVERAGE", variable="SNOW", space_agg="BY_BASIN")
conf = Config(CustomOutput=co)
```

Another similar example is the `:HRUs` command, which is rendered as a table of HRU properties. Each entry in the table is an `HRU` object. Again, HRUs can be instantiated from dictionaries or from class instances, use the style you prefer:

```{code-cell} ipython3

  conf = Config(HRUs=
  [dict(hru_id=1, area=1000, elevation=300,  latitude=50,
    longitude=-109),
  rc.HRU(hru_id=2, area=200, elevation=350,  latitude=50.5,
    longitude=-109.2)
    ])
  print(conf.rvh)
```

Gauge, StationForcing, GriddedForcing and ObservationData Commands
------------------------------------------------------------------

Meteorological forcing inputs and streamflow observations are specified using the Commands `Gauge`, `StationForcing` or `GriddedForcing`, and `ObservationData`. Ravenpy only supports reading and writing time series stored in netCDF. To facilitate the configuration of `rvt` files, each includes a class method `from_nc` that tries to infer the values of configuration options from the file's content. For example, `ObservationData.from_nc(fn, station_idx=1, alt_names=["q_obs"])` will open netCDF file `fn`, look for a variable called `q_obs` and return a `:ObservationData` command, itself holding a `ReadFromNetCDF` command with `FileNameNC`, `VarNameNC`, `DimNamesNC`, `StationIdx`, etc. Any value that should be part of the command but is not correctly extracted can be specified explicitly using extra keyword arguments given to `from_nc`.

Note that in the case of `Gauge.from_nc`, extra keyword arguments for `:Data` and `:ReadFromNetCDF` are set via the `data_kwds` dictionary, keyed by data type. The keyword `"ALL"` signifies that the keywords should be set for all variables. For example, if the forcing file does not include the gauge longitude and latitude, it could be set with `Gauge.from_nc([pr.nc, tas.nc], data_type=["PRECIP", "AVE_TEMP"], data_kwds={"ALL": {"Latitude": 45, "Longitude"=-80}})`. Linear transformations needed to convert units are inferred automatically whenever possible, if those fail, set them explicitly with `data_kwds`, e.g. `data_kwds={"PRECIP": {"Deaccumulate": True, "TimeShift": -0.25, "LinearTransform": {"scale": 1000}}}`. Note that automatic units conversion will fail for precipitation accumulated during the day which need the `:Deaccumulate` option.


Input validation
----------------

Many Raven commands can only take specific values, for example, `:PotentialMeltMethod` can take one of nine values: "POTMELT_DEGREE_DAY", "POTMELT_DATA", "POTMELT_RESTRICTED", etc. Valid options are defined in `ravenpy.config.options.PotentialMeltMethod` as an `enum.Enum` object. If `Enum` objects are not familiar, think of them as a mapping between keys and values. The full list of options for `PotentialMeltMethod` can be displayed by converting the Enum to a list:

```{code-cell} ipython3
from ravenpy.new_config import options as o

list(o.PotentialMeltMethod)
```

`Config` understands both the actual name of the option as a string, or the `Enum` option. For example, both styles below are valid:

```{code-cell} ipython3
Config(PotentialMeltMethod=o.PotentialMeltMethod.DEGREE_DAY,
       RainSnowFraction="RAINSNOW_DATA")
```
Of course, setting an option with an invalid value will raise a `ValidationError`. For example:

```{code-cell} ipython3
   Config(Routing="ROUTE_DIFFUSIVE")
```

The error message mentions what options are available to the `:Routing:` option: we should have written 'ROUTE_DIFFUSIVE_WAVE'.

This validation mechanism is used throughout the configuration, and plays dual roles: trying to convert inputs into the expected type, and it that fails, raising a `ValidationError`.


Accessing and modifying configuration options
----------------------------------------------

The configuration options are stored in the `Config` instance as attributes. To respect Python style conventions, all attributes are lower case. The rule we followed is to use underscores to split words, so that in the first example above, the `:StartDate` option is stored a `start_date`:

```{code-cell} ipython3
  print(conf.start_date)
  conf.start_date = "2023-04-01"
  print(conf.start_date)
  ```

Experience suggests that modifying configuration classes can create confusion. It is usually preferable to create and modify copies and leave the original configuration intact. This can be


`Config` supports the main configuration options, but some options might not yet be implemented. Consult section `Configuring new Commands` to learn how to add options to `ravenpy`.


Raven emulators
---------------

Raven emulators are subclasses of `Config` with various defaults set to reproduce the behavior of known hydrological models. This includes hydrological processes, the soil model, soil profiles, land-use, soil and vegetation classes, etc. Emulators packaged with ravenpy are found in `ravenpy.config.emulators`. Each one is essentially a class defining new default values for the necessary Raven commands. Attributes are given a pydantic annotation and a default using `pydantic.Field` with the alias set to the name of the Raven Command:

```{code-cell} ipython3
  class MyEmulator(Config):
    evaporation: o.Evaporation = Field(default="HARGREAVES", alias="Evaporation")
```

In the example above, the default `evaporation` attribute for `MyEmulator` is set to "HEARGREAVES". It can however be changed when instantiating the model by giving it another value, e.g. `MyEmulator(Evaporation="OUDIN")`. The annotation makes sure that whatever value is given is one of the allowed evaporation values defined in `ravenpy.config.options.Evaporation`. The alias plays two roles: 1. it allows users to define evaporation using either the Raven command name `Evaporation`, or the python attribute name `evaporation`; 2. it tells the Config that this attribute is a Raven command that should be rendered as `:Evaporation <value>`.

Note that in general, even once completely defined, emulators cannot be run as is, as they have no default parameter values, the default HRU has area, latitude and longitude set to zero, they have no meteorological forcings preset and initial conditions may be far off from reasonable values.


Symbolic expressions in emulator configuration
----------------------------------------------

`ravenpy` supports symbolic expressions in the definition of configuration parameters. That it, is is possible to define the value of Raven commands based on the value of parameters that are still undefined. This is done by

  1. defining a `dataclass` for the parameters with default values set to `pymbolic`_ `Variables`;
  2. subclassing `Config`, setting the default `params` to an instance of this dataclass;
  3. using `dataclass` attributes in symbolic expressions.

For example:

```{code-cell} ipython3
  from dataclasses import dataclass
  from pydantic import Field
  from pymbolic.primitives import Variable
  from ravenpy.new_config import Sym

  @dataclass
  class Params:
    X01 = Variable("X01")

  class MySymbolicEmulator(Config):
    params: Params = Params()
    rain_snow_transition: rc.RainSnowTransition = Field(rc.RainSnowTransition(temp=Params.X01, delta=2), alias="RainSnowTransition")
```

 This class can be instantiated, e.g. `MySymbolicEmulator()` return an instance, but this instance cannot be written to disk because parameters have not been set to numerical values. Numerical values for `params` can be set at instantiation, e.g. `MySymbolicEmulator(params=[-.5])`, or using the `set_params` method:

```{code-cell} ipython3
  sym = MySymbolicEmulator()
  num = sym.set_params([-.5])
  print(num.rvi)
```





uses `pydantic <https://docs.pydantic.dev/>`_ to
.. _`pymbolic` https://documen.tician.de/pymbolic/
