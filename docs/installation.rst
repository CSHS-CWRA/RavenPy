.. highlight:: shell

============
Installation
============


Stable release
--------------

Because RavenPy relies on the `Raven <http://raven.uwaterloo.ca>`_ binary (which is downloaded and
compiled during the setup), it must be installed in a virtual
environment (either a classic `Python venv
<https://docs.python.org/3/tutorial/venv.html>`_ or a `Conda
environment
<https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
(which must be different than the ``base`` environment).

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

For development purposes, you can add
``--install-option="--with-testdata"`` to the above ``pip install``
command, which will also download a snapshot of the the `Raven Test
Data <https://github.com/Ouranosinc/raven-testdata>`_ into your venv.


From sources
------------

The sources for RavenPy can be obtained from the GitHub repo:

.. code-block:: console

    $ git clone git://github.com/CSHS-CWRA/ravenpy

You can then follow the instructions above to setup your venv, but install with:

.. code-block:: console

   $ cd /path/to/ravenpy
   $ pip install . --verbose --install-option="--with-raven"

If you want to run the tests, install with:

.. code-block:: console

   $ cd /path/to/ravenpy
   $ pip install . --verbose --install-option="--with-raven" --install-option="--with-testdata"
   $ pip install -r requirements_dev

You can then run the test suite by doing:

.. code-block:: console

   $ pytest
