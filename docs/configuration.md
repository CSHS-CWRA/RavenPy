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

(conf)=
# Model Configuration


Raven lets hydrologists customize how models are defined, and provides examples of model configurations that can reproduce outputs from known hydrological models: HBV-EC, HMETS, Sacramento-SMA, etc. RavenPy can be used to create Python objects representing those emulators, facilitating the setup of complex modelling experiments. The following is a walk-through covering the basics of emulator configuration.

## The `Config` class

RavenPy writes configuration files using the `Config` class. `Config` recognizes most Raven commands, and has mechanisms to parse and validate the inputs, and then write them in the appropriate RV file. For example:

```{code-cell} ipython3
from ravenpy.config import Config
conf = Config(StartDate="2023-03-31")
conf.rvi
```

`Config` has recognized `StartDate` as a Raven command, and knows it should appear in the `rvi` file as a line starting with `:StartDate` followed by a date in ISO format. `StartDate` could equally have been given as a `datetime.date` or `datetime.datetime` object, and `Config` would have parsed it correctly.

Many other Raven commands are known to `Config` -- to see what commands are supported, and the type of inputs they expect, consult the class docstring with `help(Config)`. All these commands are set to a default value of `None`, and won't appear in RV files unless they are given actual values. Note that if you plan on using the `OutputReader`, make sure you set `WriteNetcdfFormat` to True.

Some Raven commands take string values, others floats, integers or dates, but some also expect lists of arguments, such as `:EvaluationMetrics`:

```{code-cell} ipython3
Config(EvaluationMetrics=["NASH_SUTCLIFFE", "RMSE"]).rvi
```

Some commands require more complex structures, for example, the configuration for  `CustomOutput` is given as a dictionary with arguments `time_per`, `stat`, `variable` and `space_agg`:

```{code-cell} ipython3
Config(CustomOutput=[{"time_per": "MONTHLY", "stat": "AVERAGE", "variable": "SNOW", "space_agg": "BY_HRU"}]).rvi
```

All Raven commands with inputs more complex than single values or simple lists are defined in `ravenpy.config.commands`. Their names match exactly with the Raven command names described in the Raven documentation. Consult the docstring to find out how each should be instantiated. Attributes can be given as dictionaries that `Config` will parse, as above, or as `Command` instances:

```{code-cell} ipython3
from ravenpy.config import commands as rc

rst = rc.RainSnowTransition(temp=-.5, delta=1)
Config(RainSnowTransition=rst).rvp
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

## Gauge, StationForcing, GriddedForcing and ObservationData commands

Meteorological forcing inputs and streamflow observations are specified using the commands `Gauge`, `StationForcing` or `GriddedForcing`, and `ObservationData`. Ravenpy only supports reading and writing time series stored in netCDF. To facilitate the configuration of `rvt` files, each includes a class method `from_nc` that tries to infer the values of configuration options from the file's content. The inference rules are based on the CF-Convention, making the configuration of models using input files that are CF compliant substantially simpler. For example, the lookup for a given forcing type (PRECIP, RAINFALL, SNOWFALL, AVE_TEMP, etc) starts with the CF standard name (`ravenpy.config.conventions.CF_RAVEN`), but if that fails, `from_nc` will try alternative names given in the `alt_names` function argument.

For example, `ObservationData.from_nc(fn, station_idx=1, alt_names=["q_obs"])` opens a netCDF file `fn`, looks for a variable named `water_volume_transport_in_river_channel`, and when that fails, looks for `q_obs`. It returns an `rc.ObservationData` instance, itself holding a `ReadFromNetCDF` command with `FileNameNC`, `VarNameNC`, `DimNamesNC`, `StationIdx`, etc. Any value that should be part of the command but is not correctly extracted can be specified explicitly using extra keyword arguments given to `from_nc`.

Note that in the case of `Gauge.from_nc`, extra keyword arguments for `:Data` and `:ReadFromNetCDF` are set via the `data_kwds` dictionary, keyed by data type. The keyword `"ALL"` signifies that the keywords should be set for all variables. For example, if the forcing file does not include the gauge longitude and latitude, it could be set with `Gauge.from_nc([pr.nc, tas.nc], data_type=["PRECIP", "AVE_TEMP"], data_kwds={"ALL": {"Latitude": 45, "Longitude"=-80}})`. Linear transformations needed to convert units are inferred automatically whenever possible. If those fail, set them explicitly with `data_kwds`, e.g. `data_kwds={"PRECIP": {"Deaccumulate": True, "TimeShift": -0.25, "LinearTransform": {"scale": 1000}}}`. Note that automatic units conversion will fail for precipitation accumulated during the day which need the `:Deaccumulate` option.


## Input validation

Many Raven commands can only take specific values, for example, `:PotentialMeltMethod` can take one of nine values: "POTMELT_DEGREE_DAY", "POTMELT_DATA", "POTMELT_RESTRICTED", etc. Valid options are defined in `ravenpy.config.options.PotentialMeltMethod` as an `enum.Enum` object. If `Enum` objects are not familiar, think of them as a static mapping between keys and values. The full list of options for `PotentialMeltMethod` can be displayed by converting the Enum to a list:

```{code-cell} ipython3
from ravenpy.config import options as o

list(o.PotentialMeltMethod)
```

`Config` understands both the option value, and the `Enum` key. For example, both styles below are valid:

```{code-cell} ipython3
Config(
       PotentialMeltMethod=o.PotentialMeltMethod.DEGREE_DAY,
       RainSnowFraction="RAINSNOW_DATA"
       ).rvi
```
Of course, setting an option with an invalid value will raise a `ValidationError`. For example, the error message below suggests there is a typo in the routing method, it should read `ROUTE_DIFFUSIVE_WAVE`:

```{code-cell} ipython3
:tags: ["raises-exception"]

Config(Routing="ROUTE_DIFFUSIVE")
```

This validation mechanism is used throughout the configuration, and will catch configuration errors well before executing the model. It relies on [pydantic](https://docs.pydantic.dev/), which compares attribute values to their type annotation. Note that types are not strict, and pydantic tries to cast the values given into the annotated type. A `ValidationError` is raised if this fails. This explains that even though `RunName` has a `str` annotation, it is still possible to pass it an integer, and it'll be converted to a string:

```{code-cell} ipython3
Config(RunName=1).rvi
```

## Accessing and modifying configuration options

The configuration options are stored in the `Config` instance as attributes. To respect Python style conventions, all attributes are lower case. The rule we followed is to use underscores to split words, so that in the first example above, the `:StartDate` option is stored as attribute `start_date`:

```{code-cell} ipython3
conf.start_date = "2023-04-01"
conf.start_date
```

While modifying configuration in place is sometimes useful, experience suggests it can create confusion and workflows that are harder to understand. We recommend instead to create modified copies of the original configuration using the `duplicate` method:

```{code-cell} ipython3
new = conf.duplicate(RunName="new", Duration=10)
print(new.rvi)
```

`duplicate` takes keyword arguments that it uses them to create a modified copy of the original configuration, which is left intact.


## Limitations

Some Raven commands are not yet supported by ravenpy. Trying to give unrecognized configuration attributes will raise a `ValidationError` saying 'extra fields not permitted'. If this happens, please submit a [feature request](https://github.com/CSHS-CWRA/RavenPy/issues).


## Raven emulators

Raven emulators are pre-configured models meant to reproduce almost exactly the behavior of known hydrological models. Emulators packaged with ravenpy are found in `ravenpy.config.emulators`. Each emulator is just a subclass of `Config`, with default values set for hydrological processes, the soil model, soil profiles, land-use, soil and vegetation classes, etc.

If you look at emulators' attributes, you'll see they are given a pydantic annotation and a default value using `pydantic.Field`, with an alias set to the name of the Raven command, for example:

```{code-cell} ipython3
from pydantic import Field

class TestEmulator(Config):
  evaporation: o.Evaporation = Field(default="PET_HARGREAVES", alias="Evaporation")
```

In the example above, the default `evaporation` attribute for `TestEmulator` is set to "PET_HEARGREAVES". It can however be changed when instantiating the model by giving it another value, e.g. `TestEmulator(Evaporation="PET_OUDIN")`. The annotation makes sure that whatever value is given is one of the allowed evaporation values defined in `ravenpy.config.options.Evaporation`.

The `alias` plays two roles:
1. it allows users to define evaporation using either the Raven command name `Evaporation`, in addition to the attribute name `evaporation`;
2. it tells the Config that this attribute is a Raven command that should be rendered as `:Evaporation <value>`.

Note that in general, even once completely defined, emulators cannot be run as is, because: (1) they have no default parameter values, (2) the default HRU has area set to zero, (3) there are no meteorological forcings, and (4) initial conditions may be far off from reasonable values.

## Symbolic expressions in emulator configuration

`ravenpy` supports symbolic expressions in the definition of configuration parameters. That is, it is possible to define the value of Raven commands based on the value of parameters that are still undefined. This is done by:

  1. defining a `dataclass` for the parameters, with type annotations allowing for symbolic variables and floats, and with default values set to [pymbolic](https://documen.tician.de/pymbolic/) `Variables`;
  2. subclassing `Config`, setting the default `params` to an instance of this dataclass;
  3. using parameters in symbolic expressions for Raven commands.

For example:

```{code-cell} ipython3
from typing import Union
from pydantic import ConfigDict
from pydantic.dataclasses import dataclass
from pymbolic.primitives import Variable
from ravenpy.config import Sym

@dataclass(config=ConfigDict(arbitrary_types_allowed=True))
class P:
  X01: Union[Variable, float] = Variable("X01")

class MyEmulator(Config):
  params: P = P()
  rain_snow_transition: rc.RainSnowTransition = Field(
  default=rc.RainSnowTransition(temp=P.X01, delta=2),
  alias="RainSnowTransition")
```

This class can be instantiated with `MyEmulator()`, but it cannot be written to disk because parameters have not been set to numerical values. Numerical values for `params` can be set at instantiation, e.g. `MyEmulator(params=[-.5])`, by attribute assignment (e.g. `conf.params=[-.5]`), or using a special-purpose `set_params` method that returns a new configuration object where symbolic expressions have been converted to numerical values.

```{code-cell} ipython3
sym = MyEmulator()
num = sym.set_params([-.5])
num
```

Such symbolic model configuration is absolutely essential for model calibration, where the calibration algorithm needs to modify parameters repeatedly. The pre-packaged emulators made available in ravenpy are already setup to perform calibration using a symbolic model configuration.
