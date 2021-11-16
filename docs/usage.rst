=====
Usage
=====

RavenPy is designed to make RAVEN easier to use from an interactive Python programming environment. It creates a temporary directory for execution, writes configuration files and creates symbolic links to input files before executing the model. Once the simulation is complete, RavenPy reads output files and converts them into native Python objects whenever possible.

Running Raven from existing configuration files
-----------------------------------------------

In a situation where Raven configuration files (`.rv*` files) exist already, Raven can be run by:

.. code-block:: python

   from ravenpy.models import Raven
   model = ravenpy.models.Raven(workdir="/tmp/testrun")
   model.configure([<list of .rv* files>])
   model([<list of forcing files>])


The `configure` method copies the configuration files to a directory, creates a symbolic link to the Raven executable in the same directory and runs it. The simulated hydrograph can then be accessed using the `q_sim` attribute, for example:

.. code-block:: python

   model.q_sim.plot()

The real value of RavenPy is however in its templating capability. Formatting tags inserted in configuration files can be replaced by Python objects injected at run time. So for example, instead of specifying the start and end dates of a simulation in the configuration file, we can write a configuration file with::

  :StartDate             {start_date}
  :EndDate               {end_date}

and then pass `start_date` and `end_date` as arguments to the `model` call:

.. code-block:: python

  import datetime as dt
  model = ravenpy.models.Raven()
  model.configure([<list of .rv* file templates>])
  model(start_date=dt.datetime(2020, 1, 1), end_date=dt.datetime(2020, 2, 1))

This templating mechanism has been put in place for all four emulated models offered by RavenPy.

Running a Raven emulated model
------------------------------

RavenPy supports pre-defined Raven-emulated models:

  - GR4J-Cemaneige
  - HMETS
  - HBV-EC
  - MOHYSE

For each one of these, `.rv` files are provided that reproduce almost perfectly the behavior of the original models and let hydrologists template typical model options. Running a simulation from an emulated model minimally requires passing a vector of model parameters and mandatory model options, as well as the list of input forcing files. RavenPy will then fill the `.rv` files with user-defined parameters and launch the simulation, for example:

.. code-block:: python

   from ravenpy.models import GR4JCN
   model = GR4JCN()
   model(ts=<list of input forcing files>,
         params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
         start_date=dt.datetime(2000, 1, 1),
         end_date=dt.datetime(2002, 1, 1),
         area=4250.6,
         elevation=843.0,
         latitude=54.4848,
         longitude=-123.3659
         )
   model.q_sim.plot()

The model configuration can be found as a zip archive in:

.. code-block:: python

   model.outputs["rv_config"]


Setting initial conditions
--------------------------
Each emulated model defines default initial conditions for its state variables (e.g. storage). Initial conditions can be set explicitly by passing the `HRUStateVariables` parameter when calling the model:

.. code-block:: python

   from ravenpy.models import GR4JCN
   from ravenpy.models.state import HRUStateVariables
   model = GR4JCN()
   model(ts=ts, hru_state=HRUStateVariables(soil0=100), **kwargs)


Resuming from a previous run
----------------------------
Once a first simulation has completed, it's possible to initialize a second simulation using the state at the end of the first simulation. This can be done from a saved `rvc` *solution* file:

.. code-block:: python

   model = GR4JCN()
   rvc = open(<path to solution.rvc>)
   model.resume(rvc)
   model(ts=ts, **kwargs)

or if a model instance already exists, simply by calling the `resume` method on it:

.. code-block:: python

   model = GR4JCN()
   model(ts=ts, start_date=dt.datetime(2000, 1, 1), end_date=dt.datetime(2002, 2, 1), **kwargs)
   model.resume()
   model(ts=ts, start_date=dt.datetime(2000, 2, 1), end_date=dt.datetime(2002, 3, 1))
