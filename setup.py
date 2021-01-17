#!/usr/bin/env python

"""The setup script."""

import shutil
import subprocess
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
    "click",
    "gdal",
    "matplotlib",
    "rasterio",
    "rioxarray",
    "statsmodels",
    "xarray",
    "xclim",
    "wheel",
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


class InstallExternalDeps(install):
    """
    Custom handler for the 'install' command, to download, extract and compile
    the source code of Raven and OSTRICH and copy the resulting binaries in a location
    available on the PATH.
    """

    user_options = install.user_options + [
        # The format is (long option, short option, description).
        ("with-binaries", None, "Download Raven and OSTRICH sources and compile them."),
    ]

    def initialize_options(self):
        """Set default values for options."""
        # Each user option must be listed here with their default value.
        install.initialize_options(self)
        self.with_binaries = False

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
            raise RuntimeError(f"There was an error while compiling {name}")

        # Copy binary binary in a location which should be available on the PATH
        # Note that self.install_scripts should correspond to <venv>/bin or ~/.local/bin,
        # depending on the install context (this is to be validated)
        target_bin_path = Path(self.install_scripts) / name
        print(f"Copying binary to {target_bin_path}")
        shutil.copy(
            self.external_deps_path / rev_name / binary_name,
            target_bin_path,
        )

    def run(self):
        if self.with_binaries:
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

        # This works with python setup.py install, but produces this error with pip install:
        # ERROR: ravenpy==0.1.0 did not indicate that it installed an .egg-info directory. Only setup.py projects generating .egg-info directories are supported.
        # super().do_egg_install()

        # This works with pip install, but has the problem that it ignores install_requires
        # when running with `python setup.py install`:
        # https://stackoverflow.com/questions/21915469/python-setuptools-install-requires-is-ignored-when-overriding-cmdclass
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
    include_package_data=True,  # This needs to be True for the files in ravenpy/data to be included
    keywords="ravenpy",
    name="ravenpy",
    packages=find_packages(
        include=[
            "ravenpy",
            "ravenpy.*",
        ],
    ),
    test_suite="tests",
    tests_require=test_requirements,
    extras_require={
        "docs": docs_requirements,
        "dev": dev_requirements,
    },
    url="https://github.com/CSHS-CWRA/ravenpy",
    version="0.1.6",
    zip_safe=False,
    cmdclass={"install": InstallExternalDeps},
)
