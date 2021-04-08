Welcome to RavenPy's documentation!
===================================


RavenPy is a Python wrapper for Raven_, accompanied by utility functions that facilitate model configuration, calibration, and evaluation.

Raven_ is an hydrological modeling framework that lets hydrologists build hydrological models by combining different hydrological processes together. It can also be used to emulate a variety of existing lumped and distributed models. Model structure, parameters, initial conditions and forcing files are configured in text files, which Raven parses to build and run hydrological simulations. A detailed description about modeling capability of Raven_ can be found in the `docs`_.

RavenPy provides a Python interface to Raven_, automating the creation of configuration files and allowing the model to be launched from Python. Results, or errors, are automatically parsed and exposed within the programming environment. This facilitates the launch of parallel simulations, multi-model prediction ensembles, sensitivity analyses and other experiments involving a large number of model runs.

The code is currently undergoing a number of changes to support semi-distributed watersheds, so expect some API changes over the next versions.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   usage
   notebooks/index
   contributing
   authors
   history

.. toctree::
   :maxdepth: 1
   :caption: Utility Scripts

   scripts

.. toctree::
   :maxdepth: 1
   :caption: User API

   user_api

.. toctree::
   :maxdepth: 1
   :caption: All Modules

   modules


Credits
=======
RavenPy's development has been funded by CANARIE_ and Ouranos_ and would be not be possible without the help of Juliane Mai and James Craig.

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _`Raven`: http://raven.uwaterloo.ca
.. _`CANARIE`: https://www.canarie.ca
.. _`docs`: http://raven.uwaterloo.ca/files/v3.0.1/RavenManual_v3.0.1.pdf
.. _`Ouranos`: http://ouranos.ca
