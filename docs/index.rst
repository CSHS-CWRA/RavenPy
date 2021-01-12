Welcome to RavenPy's documentation!
===================================

RavenPy is a Python wrapper around the Raven_ hydrological modeling framework. Raven lets hydrologists build hydrological models by combining different hydrological processes together. It can also be used to emulate a variety of existing lumped and distributed models. Model structure, parameters, initial conditions and forcing files are configured in text files, which Raven parses to build and run hydrological simulations. A detailed description about modeling capability and configuration options can be found in the `docs`_.

RavenPy downloads and compiles the most recent version of the Raven, and offers utilities that can automatically write model configuration files, run the executable and parse output files. This let's scientists work with Raven within a Python programming environment, facilitating the launch of parallel simulations, multi-model prediction ensembles, sensitivity analyses and other experiments involving a large number of model runs.

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
   :caption: Scripts

   scripts


.. toctree::
   :maxdepth: 1
   :caption: API

   user_api


Credits
=======
RavenPy's development has been funded by CANARIE_.

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _`Raven`: http://raven.uwaterloo.ca
.. _`CANARIE`: https://www.canarie.ca
.. _`docs`: http://raven.uwaterloo.ca/files/v3.0.1/RavenManual_v3.0.1.pdf
