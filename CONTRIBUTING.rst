.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/CSHS-CWRA/RavenPy/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement" and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

RavenPy could always use more documentation, whether as part of the official RavenPy docs, in docstrings, or even on the web in blog posts, articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/CSHS-CWRA/RavenPy/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome. :)

Get Started!
------------

.. note::

    If you are new to using `GitHub <https://github.com/>`_ and ``git``, please read `this guide <https://guides.github.com/activities/hello-world/>`_ first.

.. warning::

    Anaconda Python users: Due to the complexity of some packages, the default dependency solver can take a long time to resolve the environment. Consider running the following commands in order to speed up the process:

        .. code-block:: console

            conda install -n base conda-libmamba-solver
            conda config --set solver libmamba

        For more information, please see the following link: https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community

    Alternatively, you can use the `mamba <https://mamba.readthedocs.io/en/latest/index.html>`_ package manager, which is a drop-in replacement for ``conda``. If you are already using `mamba`, replace the following commands with ``mamba`` instead of ``conda``.

Ready to contribute? Here's how to set up `ravenpy` for local development.

#. First, clone the ``RavenPy`` repo locally.

    * If you are not a ``RavenPy`` collaborator, first fork the ``RavenPy`` repo on GitHub, then clone your fork locally.

        .. code-block:: console

            git clone git@github.com:your_name_here/RavenPy.git

    * If you are a ``RavenPy`` collaborator, clone the ``RavenPy`` repo directly.

        .. code-block:: console

            git clone git@github.com:CSHS-CWRA/RavenPy.git

#. Install your local copy into a development environment. You can create a new Anaconda development environment with:

    .. code-block:: console

        conda env create -f environment-dev.yml
        conda activate ravenpy-dev
        (ravenpy-dev) make dev

    If you are on Windows, replace the ``make dev`` command with the following:

        .. code-block:: console

            python -m pip install -e .[dev]
            pre-commit install

    This installs ``ravenpy`` in an "editable" state, meaning that changes to the code are immediately seen by the environment. To ensure a consistent coding style, `make dev` also installs the ``pre-commit`` hooks to your local clone.

    On commit, ``pre-commit`` will check that ``black``, ``blackdoc``, ``isort``, ``flake8``, and ``ruff`` checks are passing, perform automatic fixes if possible, and warn of violations that require intervention. If your commit fails the checks initially, simply fix the errors, re-add the files, and re-commit.

    You can also run the hooks manually with:

        .. code-block:: console

            pre-commit run -a

    If you want to skip the ``pre-commit`` hooks temporarily, you can pass the `--no-verify` flag to `git commit`.

#. Create a branch for local development:

    .. code-block:: console

        git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

+#. When you're done making changes, we **strongly** suggest running the tests in your environment or with the help of ``tox``:


    .. code-block:: console

        make lint
        python -m pytest
        # Or, to run multiple build tests
        python -m tox

#. Commit your changes and push your branch to GitHub:

    .. code-block:: console

        git add .
        git commit -m "Your detailed description of your changes."
        git push origin name-of-your-bugfix-or-feature

    If ``pre-commit`` hooks fail, try fixing the issues, re-staging the files to be committed, and re-committing your changes (or, if need be, you can skip them with `git commit --no-verify`).

#. Submit a `Pull Request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_ through the GitHub website.

#. When pushing your changes to your branch on GitHub, the documentation will automatically be tested to reflect the changes in your Pull Request. This build process can take several minutes at times. If you are actively making changes that affect the documentation and wish to save time, you can compile and test your changes beforehand locally with:

    .. code-block:: console

        # To generate the html and open it in your browser
        make docs
        # To only generate the html
        make autodoc
        make -C docs html
        # To simply test that the docs pass build checks
        python -m tox -e docs

#. If changes to your branch are made on GitHub, you can update your local branch with:

    .. code-block:: console

        git checkout name-of-your-bugfix-or-feature
        git fetch
        git pull origin name-of-your-bugfix-or-feature

    If you have merge conflicts, you might need to replace `git pull` with `git merge` and resolve the conflicts manually.
    Resolving conflicts from the command line can be tricky. If you are not comfortable with this, you can ignore the last command and instead use a GUI like PyCharm or Visual Studio Code to merge the remote changes and resolve the conflicts.

#. Before merging, your Pull Request will need to be based on the `main` branch of the ``RavenPy`` repository. If your branch is not up-to-date with the `main` branch, you can perform similar steps as above to update your branch:

    .. code-block:: console

        git checkout name-of-your-bugfix-or-feature
        git fetch
        git pull origin main

    See the previous step for more information on resolving conflicts.

#. Once your Pull Request has been accepted and merged to the `main` branch, several automated workflows will be triggered:

    - The ``bump-version.yml`` workflow will automatically bump the patch version when pull requests are pushed to the `main` branch on GitHub. **It is not recommended to manually bump the version in your branch when merging (non-release) pull requests (this will cause the version to be bumped twice).**
    - `ReadTheDocs` will automatically build the documentation and publish it to the `latest` branch of `ravenpy` documentation website.
    - If your branch is not a fork (i.e. you are a maintainer), your branch will be automatically deleted.

You will have contributed to ``ravenpy``!

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

#. The pull request should include tests and should aim to provide `code coverage <https://en.wikipedia.org/wiki/Code_coverage>`_ for all new lines of code. You can use the `--cov-report html --cov ravenpy` flags during the call to ``pytest`` to generate an HTML report and analyse the current test coverage.

#. All functions should be documented with `docstrings` following the `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html>`_ format.

#. If the pull request adds functionality, either update the documentation or create a new notebook that demonstrates the feature. Library-defining features should also be listed in ``README.rst``.

#. The pull request should work for all currently supported Python versions. Check the `pyproject.toml` or `tox.ini` files for the list of supported versions.

Tips
----

To run a subset of tests:

.. code-block:: console

    python -m pytest tests/test_ravenpy.py

You can also directly call a specific test class or test function using:

.. code-block:: console

    python -m pytest tests/test_ravenpy.py::TestClassName::test_function_name

For more information on running tests, see the `pytest documentation <https://docs.pytest.org/en/latest/usage.html>`_.

To run specific code style checks:

.. code-block:: console

    python -m black --check src/ravenpy tests
    python -m isort --check src/ravenpy tests
    python -m blackdoc --check src/ravenpy docs
    python -m ruff check src/ravenpy tests
    python -m flake8 src/ravenpy tests
    validate-docstrings src/ravenpy/**.py

To get ``black``, ``isort``, ``blackdoc``, ``ruff``, ``flake8`` (with the ``flake8-rst-docstrings`` plugin), and ``numpydoc`` (for ``validate-docstrings``), simply install them with ``pip`` (or ``conda``) into your environment.

Translations
------------

If you would like to contribute to the French translation of the documentation, you can do so by running the following command:

    .. code-block:: console

        make initialize-translations

This will create or update the French translation files in the `docs/locales/fr/LC_MESSAGES` directory. You can then edit the `.po` files in this directory to provide translations for the documentation.

Code of Conduct
---------------

Please note that this project is released with a `Contributor Code of Conduct <https://github.com/CSHS-CWRA/RavenPy/blob/master/CODE_OF_CONDUCT.md>`_.
By participating in this project you agree to abide by its terms.
