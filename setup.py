#!/usr/bin/env python

"""The setup script."""

import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path
from urllib.parse import urljoin

from setuptools import Distribution, find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install

RAVEN_ZIP_NAME = "Raven-rev318"
OSTRICH_GIT_VERSION = "21.03.16"

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "click",
    "climpred>=2.1",
    "dask",
    "haversine",
    "matplotlib",
    "netCDF4",
    "numpy",
    "pandas",
    "pydantic",
    "requests",
    "scipy",
    "setuptools<58.0.0",
    "statsmodels",
    "xarray>=0.18",
    "cf-xarray",
    "xclim>=0.26.1",
    "wheel",
    "xskillscore",
]

test_requirements = [
    "pytest>=3",
]

docs_requirements = [
    dependency for dependency in open("requirements_docs.txt").readlines()
]

gis_requirements = [
    dependency for dependency in open("requirements_gis.txt").readlines()
]
# Special GDAL handling
try:
    gdal_version = subprocess.run(
        ["gdal-config", "--version"], capture_output=True
    ).stdout.decode("utf-8")
    gis_requirements.append(f"gdal=={gdal_version}")
except subprocess.CalledProcessError:
    pass

dev_requirements = gis_requirements.copy()
dev_requirements.extend(
    [dependency for dependency in open("requirements_dev.txt").readlines()]
)


# Idea taken from: https://stackoverflow.com/a/25176606/787842
class OnlyGetScriptPath(install):
    def run(self):
        # does not call install.run() by design
        self.distribution.install_scripts = self.install_scripts  # noqa


def get_setuptools_install_scripts_dir():
    dist = Distribution({"cmdclass": {"install": OnlyGetScriptPath}})
    dist.dry_run = True  # not sure if necessary, but to be safe
    dist.parse_config_files()
    command = dist.get_command_obj("install")
    command.ensure_finalized()
    command.run()
    return dist.install_scripts


def create_external_deps_install_class(command_cls):
    """
    Class factory command to implement the customized binary download + compile + install logic
    for both the install and develop command contexts.
    """

    class InstallExternalDeps(command_cls):
        """
        Custom handler for the 'install' and 'develop' commands, to download, extract and compile
        the source code of Raven and OSTRICH and copy the resulting binaries in a location
        available on the PATH.
        """

        user_options = command_cls.user_options + [
            # The format is (long option, short option, description).
            (
                "with-binaries",
                None,
                "Download Raven and OSTRICH sources and compile them.",
            ),
        ]

        def initialize_options(self):
            """Set default values for options."""
            # Each user option must be listed here with their default value.
            command_cls.initialize_options(self)
            self.with_binaries = False  # noqa

        def finalize_options(self):
            command_cls.finalize_options(self)

        def install_binary_dep(
            self, url, name, rev_name, binary_name, make_target="", src_folder=None
        ):
            print(f"Downloading {name} source code..")
            print(
                f"{urljoin(url, rev_name)}.zip", self.external_deps_path / f"{name}.zip"
            )

            urllib.request.urlretrieve(
                f"{urljoin(url, rev_name)}.zip", self.external_deps_path / f"{name}.zip"
            )

            print(f"Extracting {name} source code..")
            with zipfile.ZipFile(
                self.external_deps_path / f"{name}.zip", "r"
            ) as zip_ref:
                zip_ref.extractall(self.external_deps_path)

            print(f"Compiling {name}..")

            folder = src_folder if src_folder else rev_name

            try:
                print(self.external_deps_path / folder)

                subprocess.check_call(
                    f"make {make_target}",
                    cwd=self.external_deps_path / folder,
                    shell=True,
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"There was an error while compiling {name}") from e

            # Copy binary in a location which should be available on the PATH
            # Note 1: if command_cls==install, self.install_scripts should correspond to <venv>/bin or ~/.local/bin
            # Note 2: if command_cls==develop, self.install_scripts is None so we are using a trick to get the value
            #         it would have with the install command
            scripts_dir = self.install_scripts or get_setuptools_install_scripts_dir()
            target_bin_path = Path(scripts_dir) / name

            print(
                f"Copying binary from "
                f"{self.external_deps_path.joinpath(src_folder if src_folder else rev_name).joinpath(binary_name)}"
            )
            print(f"To {target_bin_path}")

            shutil.copy(
                self.external_deps_path.joinpath(
                    src_folder if src_folder else rev_name
                ).joinpath(binary_name),
                target_bin_path,
            )

        def run(self):

            if self.with_binaries:
                self.external_deps_path = Path("./external_deps")  # noqa
                self.external_deps_path.mkdir(exist_ok=True)

                url = "https://www.civil.uwaterloo.ca/jmai/raven/"
                self.install_binary_dep(url, "raven", RAVEN_ZIP_NAME, "Raven.exe")

                url = "https://github.com/usbr/ostrich/archive/refs/tags/"
                self.install_binary_dep(
                    url,
                    "ostrich",
                    f"v{OSTRICH_GIT_VERSION}",
                    "Ostrich",
                    "GCC",
                    src_folder=f"ostrich-{OSTRICH_GIT_VERSION}/make",
                )

            # This works with python setup.py install, but produces this error with pip install:
            # ERROR: ravenpy==0.1.0 did not indicate that it installed an .egg-info directory.
            # Only setup.py projects generating .egg-info directories are supported.
            # super().do_egg_install()

            # This works with pip install, but has the problem that it ignores install_requires
            # when running with `python setup.py install`:
            # https://stackoverflow.com/questions/21915469/python-setuptools-install-requires-is-ignored-when-overriding-cmdclass
            command_cls.run(self)

    return InstallExternalDeps


setup(
    author="David Huard",
    author_email="huard.david@ouranos.ca",
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description="A Python wrapper to setup and run the hydrologic modelling framework Raven.",
    entry_points={
        "console_scripts": [
            "ravenpy=ravenpy.cli:main",
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
    test_suite="tests",
    tests_require=test_requirements,
    extras_require=dict(
        dev=dev_requirements,
        docs=docs_requirements,
        gis=gis_requirements,
    ),
    url="https://github.com/CSHS-CWRA/ravenpy",
    version="0.7.6",
    zip_safe=False,
    cmdclass={
        "install": create_external_deps_install_class(install),
        "develop": create_external_deps_install_class(develop),
    },
)
