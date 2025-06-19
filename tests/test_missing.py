import logging
import os
import shutil
import sys

import pytest


def filter_raven():
    # Find the real conda-installed tool path
    real_tool = shutil.which("raven")
    assert real_tool is not None, "raven must be installed for this test"

    # Get the directory the conda tool is in
    conda_tool_dir = os.path.dirname(real_tool)

    # Create a new PATH that includes everything *except* the conda env's tool directory
    filtered_path = os.pathsep.join(
        [
            p
            for p in os.environ["PATH"].split(os.pathsep)
            if os.path.abspath(p) != os.path.abspath(conda_tool_dir)
        ]
    )

    return filtered_path


@pytest.fixture
def hide_module(monkeypatch):
    def _hide(name):
        monkeypatch.setitem(sys.modules, name, None)

    return _hide


class TestMissing:

    def test_missing_raven_binary(self, monkeypatch, tmpdir, hide_module):
        """Test for behaviour when binary is missing from the system path."""

        # Set up a temporary directory to simulate the absence of the raven binary
        filtered_path = filter_raven()
        monkeypatch.setenv("PATH", f"{tmpdir}{os.pathsep}{filtered_path}")

        # Hide the raven_hydro module
        hide_module("raven_hydro")
        hide_module("raven_hydro._version")
        hide_module("raven_hydro.libraven")

        # Force the raven modules to be reloaded
        del sys.modules["ravenpy"]
        del sys.modules["ravenpy._raven"]
        del sys.modules["ravenpy.config.defaults"]

        # Now the tool should be "missing"
        # if running from tox in a conda environment, the raven binary is likely to be present
        if os.environ.get("CONDA_PREFIX") and os.environ.get("TOX"):
            logging.info(
                "Running in a development conda environment with tox. Raven is expected to be present."
            )
        else:
            assert shutil.which("raven") is None

        # Loading the module should raise a RuntimeError
        with pytest.raises(RuntimeError):
            import ravenpy  # noqa: F401

        # Check that setting the RAVENPY_RAVEN_BINARY_PATH environment variable works
        monkeypatch.setenv("RAVENPY_RAVEN_BINARY_PATH", f"some/path/to/raven")
        import ravenpy

        assert ravenpy.RAVEN_EXEC_PATH == "some/path/to/raven"

        # Check that the raven_hydro library is not imported
        from ravenpy.config.defaults import __raven_version__

        assert __raven_version__ == "0.0.0"
