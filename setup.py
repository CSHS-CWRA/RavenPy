#!/usr/bin/env python

"""The setup script."""

import os
import shutil
import subprocess
import urllib.request
import zipfile
from pathlib import Path
from typing import Optional, Union
from urllib.parse import urljoin

# Note: setuptools < 65.6 is needed for some dependencies (see: https://github.com/pypa/setuptools/issues/3693)
from setuptools import Distribution, find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install

RAVEN_VERSION = "3.6"
OSTRICH_GIT_VERSION = "21.03.16"

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "cf-xarray",
    "click",
    "climpred>=2.1",
    "dask",
    "haversine",
    "matplotlib",
    "netCDF4",
    "numpy",
    "pandas",
    "pint>=0.20",
    "pydantic",
    "requests",
    "scipy",
    "spotpy",
    "statsmodels",
    "wheel",
    "xarray<2022.11.0",  # Pinned due to incompatibility with climpred @ 2.2.0
    "xclim>=0.40.0",
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
on_conda = os.getenv("CONDA_BUILD")
if on_conda == "1":
    gis_requirements.append("gdal")
else:
    try:
        gdal_version = subprocess.run(
            ["gdal-config", "--version"], capture_output=True
        ).stdout.decode("utf-8")
        gis_requirements.append(f"gdal=={gdal_version}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        pass

dev_requirements = gis_requirements.copy()
dev_requirements.extend(
    [dependency for dependency in open("requirements_dev.txt").readlines()]
)


# Idea taken from: https://stackoverflow.com/a/25176606/787842
class OnlyGetScriptPath(install):
    def run(self):
        # does not call install.run() by design
        self.distribution.install_scripts = self.install_scripts


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

        external_deps_path = None
        with_binaries = False

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
            self.with_binaries = False

        def finalize_options(self):
            command_cls.finalize_options(self)

        def install_binary_dep(
            self,
            url,
            name: str,
            version: str,
            rev_name: Optional[str] = None,
            binary_name: str = "",
            make_target: str = "",
            src_folder: Optional[Union[str, os.PathLike]] = None,
            remove_line: Optional[str] = None,
        ):
            print(f"Downloading {name} source code..")
            if rev_name:
                file_path = f"v{(Path(version) / rev_name).as_posix()}"
            else:
                file_path = f"v{version}"

            print(
                f"{urljoin(url, file_path)}.zip",
                self.external_deps_path / f"{name}.zip",
            )
            urllib.request.urlretrieve(
                f"{urljoin(url, file_path)}.zip",
                self.external_deps_path / f"{name}.zip",
            )

            print(f"Extracting {name} source code..")
            if rev_name:
                out_folder = self.external_deps_path.joinpath(rev_name)
            else:
                out_folder = self.external_deps_path
            with zipfile.ZipFile(
                self.external_deps_path / f"{name}.zip", "r"
            ) as zip_ref:
                zip_ref.extractall(out_folder)

            print(f"Compiling {name}..")
            src_folder = src_folder if src_folder else rev_name
            c_filepath = self.external_deps_path / src_folder
            try:
                print(c_filepath)

                # Hacky patch fix until we can safely remove all this logic
                if remove_line:
                    print("Patching Makefile..")
                    with open(c_filepath.joinpath("Makefile"), "r+") as f:
                        d = f.readlines()
                        f.seek(0)
                        for i in d:
                            if remove_line not in i:
                                f.write(i)
                        f.truncate()

                subprocess.check_call(
                    f"make {make_target}",
                    cwd=c_filepath,
                    shell=True,
                )
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"There was an error while compiling {name}") from e

            # Copy binary in a location which should be available on the PATH
            # Note 1: if command_cls==install, self.install_scripts should correspond to <venv>/bin or ~/.local/bin
            # Note 2: if command_cls==develop, self.install_scripts is None, so we are using a trick to get the value
            #         it would have with the `install` command
            scripts_dir = self.install_scripts or get_setuptools_install_scripts_dir()
            target_bin_path = Path(scripts_dir) / name

            print(
                f"Copying binary from: "
                f"{self.external_deps_path.joinpath(src_folder).joinpath(binary_name)}\n"
                f"To: {target_bin_path}"
            )
            shutil.copy(
                self.external_deps_path.joinpath(src_folder).joinpath(binary_name),
                target_bin_path,
            )

        def run(self):
            if self.with_binaries:
                self.external_deps_path = Path().cwd().joinpath("external_deps")
                self.external_deps_path.mkdir(exist_ok=True)

                url = "https://www.civil.uwaterloo.ca/raven/files/"
                self.install_binary_dep(
                    url,
                    "raven",
                    version=RAVEN_VERSION,
                    rev_name=f"RavenSource_v{RAVEN_VERSION}",
                    binary_name="Raven.exe",
                    remove_line="CXXFLAGS += -c++11",
                )

                url = "https://github.com/usbr/ostrich/archive/refs/tags/"
                self.install_binary_dep(
                    url,
                    "ostrich",
                    version=OSTRICH_GIT_VERSION,
                    binary_name="Ostrich",
                    make_target="GCC",
                    src_folder=Path(f"ostrich-{OSTRICH_GIT_VERSION}/make"),
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
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
        "Topic :: Scientific/Engineering :: GIS",
        "Topic :: Scientific/Engineering :: Hydrology",
    ],
    description="A Python wrapper to setup and run the hydrologic modelling framework Raven.",
    entry_points={
        "console_scripts": ["ravenpy=ravenpy.cli:main"],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    package_data={"ravenpy": ["*.csv", "*.zip"]},
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
    version="0.11.0",
    zip_safe=False,
    cmdclass={
        "install": create_external_deps_install_class(install),
        "develop": create_external_deps_install_class(develop),
    },
)
