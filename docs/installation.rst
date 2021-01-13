.. highlight:: shell

============
Installation
============


Stable release
--------------

Because RavenPy relies on the `Raven <http://raven.uwaterloo.ca>`_ binary (which is downloaded and
compiled during the setup), it is preferable to install it in a virtual
environment (either a classic `Python venv
<https://docs.python.org/3/tutorial/venv.html>`_ or a `Conda
environment
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

As the ``GDAL`` Python library requires some system-level dependencies
to be properly compiled, you will probably have to do:

.. code-block:: console

   $ sudo apt update && sudo apt install libgdal-dev

Once your virtual has been activated, it will not hurt to do:

.. code-block:: console

   $ pip install -U pip setuptools wheel

Next you can do:

.. code-block:: console

   $ pip install ravenpy --verbose --install-option="--with-raven"

If successful, this should install a ``raven`` binary in the ``bin``
folder of your venv, which should then be already available in your
path, and that you can verify by doing:

.. code-block:: console

   $ which raven

If you don't use the ``--install-option="--with-raven"`` option or if
there's a problem with it, you can also manage the necessary binaries
yourself (see the next section).


From sources (for development and testing)
------------------------------------------

The sources for RavenPy can be obtained from the GitHub repo:

.. code-block:: console

    $ git clone git://github.com/CSHS-CWRA/ravenpy

You must download and install `Raven <http://raven.uwaterloo.ca>`_ and
`OSTRICH <http://www.civil.uwaterloo.ca/envmodelling/Ostrich.html>`_
and they must be available on your PATH. Alternatively, you can supply
``RAVENPY_RAVEN_BINARY_PATH`` and ``RAVENPY_OSTRICH_BINARY_PATH`` env
variables to your Python interpreter.

Once you have created and activated your venv, you can do:

.. code-block:: console

   $ cd /path/to/ravenpy
   $ pip install --verbose --editable .

If you want to run the tests, install the requirements:

.. code-block:: console

   $ pip install -r requirements_dev

Then clone the Raven Test Data repo somewhere on your disk:

.. code-block:: console

    $ git clone git@github.com:Ouranosinc/raven-testdata.git

You can then run the test suite by doing:

.. code-block:: console

   $ RAVENPY_TESTDATA_PATH=/path/to/raven-testdata pytest
