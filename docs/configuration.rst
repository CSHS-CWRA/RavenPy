===================
Model configuration
===================
Raven lets hydrologist customize how models are defined, and provides examples of model configuration that can reproduce outputs from known hydrological models: HBV-EC, HMETS, Sacramento-SMA, etc. RavenPy can be used to create Python objects representing those emulators, facilitating the setup of complex modeling experiments. The following describes the basics of emulator configuration.

Configuring Raven emulators
---------------------------
The real value of RavenPy lies in its capability to write configuration files from a Python environment. This is done by the `Config` class, which is able to recognize most Raven commands, and write them in the appropriate RV file. For example:

.. code-block:: python

  from ravenpy.config import Config

  conf = Config(StartDate="2023-03-31")
  print(conf.rvi)

  :StartDate            2023-03-31 00:00:00
  :WriteNetcdfFormat

The `Config` has recognized `StartDate` as a Raven command, and converted it to a line within the `rvi` file. `StartDate` could equally have been given as a `datetime.date` or `datetime.datetime` object. By default, the `Config` class sets `WriteNetcdfFormat` to True to make sure outputs can be parsed by `OutputReader`. The list of all available configuration options can be consulted with `help(Config)`.

Many Raven commands can only take specific values, for example, `:PotentialMeltMethod`. Valid options are listed in `ravenpy.config.options` as `enum.Enum` objects. The `Config` class will understand both the actual name of the option as a string, or the `Enum` option. For example, both styles below are valid:

.. code-block:: python

  from ravenpy.config import options as o
  Config(PotentialMeltMethod=o.DEGREE_DAY,
         RainSnowFraction="RAINSNOW_DATA")

What is particularly valuable is that `Config` uses `pydantic <https://docs.pydantic.dev/>`_ to compare values against expectations. If an option has a typo, it will raise a `ValidationError`.For example:

.. code-block:: python

   Config(Routing="ROUTE_DIFFUSIVE")

  ValidationError: 1 validation error for Config
  Routing
    value is not a valid enumeration member; permitted: 'ROUTE_DIFFUSIVE_WAVE', 'ROUTE_HYDROLOGIC', 'ROUTE_NONE', 'ROUTE_STORAGE_COEFF', 'ROUTE_PLUG_FLOW', 'MUSKINGUM' (type=type_error.enum; enum_values=[<Routing.DIFFUSIVE_WAVE: 'ROUTE_DIFFUSIVE_WAVE'>, <Routing.HYDROLOGIC: 'ROUTE_HYDROLOGIC'>, <Routing.NONE: 'ROUTE_NONE'>, <Routing.STORAGE_COEFF: 'ROUTE_STORAGE_COEFF'>, <Routing.PLUG_FLOW: 'ROUTE_PLUG_FLOW'>, <Routing.MUSKINGUM: 'MUSKINGUM'>])

We can rapidly realize that we should have written 'ROUTE_DIFFUSIVE_WAVE'.

The configuration options are stored in the `Config` instance as class attributes. To respect Python style conventions, all attributes are lower case. The rule we followed is to use underscores to split words, so that in the first example above, to retrieve the start date we would write:

.. code-block:: python

  conf.start_date
  cftime.DatetimeProlepticGregorian(2023, 3, 31, 0, 0, 0, 0, has_year_zero=True)
