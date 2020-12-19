#!/usr/bin/env python

"""The setup script."""

import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

from setuptools import find_packages, setup
from setuptools.command.install import install

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "Click~=7.0",
    "gdal==3.0.4",
    "rasterio",
    "statsmodels",
    "xarray==0.16",
    "xclim==0.18",
    "wheel",
]

setup_requirements = [
    "pytest-runner",
]

test_requirements = [
    "pytest>=3",
]

docs_requirements = [
    dependency for dependency in open("requirements_docs.txt").readlines()
]

dev_requirements = [
    dependency for dependency in open("requirements_dev.txt").readlines()
]


class InstallBinaryDeps(install):
    """
    Custom handler for the 'install' command, to download, extract and compile
    the source code of Raven and OSTRICH and copy the resulting binaries in the
    "bin" folder of the current venv.
    """

    user_options = install.user_options + [
        # The format is (long option, short option, description).
        ("with-raven", None, "Download Raven and OSTRICH sources and compile them."),
        ("with-testdata", None, "Download RavenPy test data used for unit tests."),
    ]

    def initialize_options(self):
        """Set default values for options."""
        # Each user option must be listed here with their default value.
        install.initialize_options(self)
        self.with_raven = False
        self.with_testdata = False

    def finalize_options(self):
        install.finalize_options(self)

    def install_binary_dep(self, url, name, rev_name, binary_name, make_target=""):
        print(f"Downloading {name} source code..")
        urllib.request.urlretrieve(
            f"{url}/{rev_name}.zip", self.external_deps_path / f"{name}.zip"
        )

        print(f"Extracting {name} source code..")
        with zipfile.ZipFile(self.external_deps_path / f"{name}.zip", "r") as zip_ref:
            zip_ref.extractall(self.external_deps_path)

        print(f"Compiling {name}..")
        try:
            subprocess.check_call(
                f"make {make_target}",
                cwd=self.external_deps_path / rev_name,
                shell=True,
            )
        except subprocess.CalledProcessError as e:
            print(e)
            exit(f"There was an error while compiling {name}")

        #  Copy binary into venv bin folder (so it should be in the path when the venv is active)
        shutil.copy(
            self.external_deps_path / rev_name / binary_name,
            Path(sys.prefix) / "bin" / name,
        )

    def run(self):
        if sys.base_prefix == sys.prefix and not os.getenv("CONDA_PREFIX"):
            exit("Error: Please install RavenPy in a virtual environment!")

        if self.with_raven:
            self.external_deps_path = Path("./external_deps")
            self.external_deps_path.mkdir(exist_ok=True)

            url = "http://www.civil.uwaterloo.ca/jmai/raven/"
            self.install_binary_dep(url, "raven", "Raven-rev288", "Raven.exe")
            self.install_binary_dep(
                url,
                "ostrich",
                "Ostrich_2017-12-19_plus_progressJSON",
                f"OstrichGCC",
                "GCC",
            )

        # this works with python setup.py install
        # super().do_egg_install()

        # this works with pip install:
        install.run(self)


setup(
    author="David Huard",
    author_email="huard.david@ouranos.ca",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="A Python wrapper to setup and run the hydrologic modelling framework Raven.",
    entry_points={
        "console_scripts": [
            "ravenpy=ravenpy.cli:cli",
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    keywords="ravenpy",
    name="ravenpy",
    packages=find_packages(
        include=[
            "ravenpy",
            "ravenpy.*",
        ],
    ),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    extras_require={
        "docs": docs_requirements,
        "dev": dev_requirements,
    },
    url="https://github.com/CSHS-CWRA/ravenpy",
    version="0.1.0",
    zip_safe=False,
    cmdclass={"install": InstallBinaryDeps},
)
