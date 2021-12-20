.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/CSHS-CWRA/ravenpy/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

RavenPy could always use more documentation, whether as part of the
official RavenPy docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/CSHS-CWRA/ravenpy/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `ravenpy` for local development.

1. Fork the `ravenpy` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/ravenpy.git

3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

    $ mkvirtualenv ravenpy
    $ cd ravenpy/
    $ python setup.py develop

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass flake8 and the
   tests, including testing other Python versions with tox::

    $ flake8 ravenpy tests
    $ black --check ravenpy tests
    $ python setup.py test or pytest
    $ tox

   To get flake8 and tox, just pip install them into your virtualenv.

6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.7, 3.8, and 3.9. Check
   https://travis-ci.com/CSHS-CWRA/ravenpy/pull_requests
   and make sure that the tests pass for all supported Python versions.

Tips
----

To run a subset of tests::

$ pytest tests.test_ravenpy


Versioning/Tagging
------------------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

$ bumpversion patch # possible: major / minor / patch
$ git push
$ git push --tags

Packaging
---------

When a new version has been minted (features have been successfully integrated test coverage and stability is adequate),
maintainers should update the pip-installable package (wheel and source release) on PyPI as well as the binary on conda-forge.

The simple approach
~~~~~~~~~~~~~~~~~~~

The simplest approach to packaging for general support (pip wheels) requires the following packages installed:
 * setuptools
 * wheel
 * twine

From the command line on your Linux distribution, simply run the following from the clone's main dev branch::

    # To build the packages (sources and wheel)
    $ python setup.py sdist bdist_wheel

    # To upload to PyPI
    $ twine upload dist/*

The new version based off of the version checked out will now be available via `pip` (`$ pip install ravenpy`).

Releasing on conda-forge
~~~~~~~~~~~~~~~~~~~~~~~~

Initial Release
^^^^^^^^^^^^^^^

In order to prepare an initial release on conda-forge, we *strongly* suggest consulting the following links:
 * https://conda-forge.org/docs/maintainer/adding_pkgs.html
 * https://github.com/conda-forge/staged-recipes

Before updating the main conda-forge recipe, we echo the conda-forge documentation and *strongly* suggest performing the following checks:
 * Ensure that dependencies and dependency versions correspond with those of the tagged version, with open or pinned versions for the `host` requirements.
 * If possible, configure tests within the conda-forge build CI (e.g. `imports: ravenpy`, `commands: pytest ravenpy`)

Subsequent releases
^^^^^^^^^^^^^^^^^^^

If the conda-forge feedstock recipe is built from PyPI, then when a new release is published on PyPI, `regro-cf-autotick-bot` will open Pull Requests automatically on the conda-forge feedstock.
It is up to the conda-forge feedstock maintainers to verify that the package is building properly before merging the Pull Request to the main branch.
