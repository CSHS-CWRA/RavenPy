import numpy as np
import pytest
import xarray as xr
from climpred import HindcastEnsemble

from ravenpy import OutputReader
from ravenpy.utilities.new_config.forecasting import (
    climatology_esp,
    hindcast_climatology_esp,
    to_climpred_hindcast_ensemble,
    warm_up,
)


def test_warm_up(minimal_emulator, tmp_path):
    conf = warm_up(minimal_emulator, duration=10, workdir=tmp_path)
    wup = OutputReader(path=tmp_path / "output")

    assert conf.start_date == minimal_emulator.start_date
    assert len(wup.hydrograph.q_sim.time) == 11
    assert conf.hru_states != minimal_emulator.hru_states


def test_climatology_esp(minimal_emulator, tmp_path):
    esp = climatology_esp(minimal_emulator, workdir=tmp_path, years=[1955, 1956])
    np.testing.assert_array_equal(esp.storage.member, [1955, 1956])
    assert len(esp.hydrograph.time) == minimal_emulator.duration + 1


def test_hindcast_climatology_esp(minimal_emulator, tmp_path, salmon_meteo):
    hc = hindcast_climatology_esp(
        minimal_emulator,
        warm_up_duration=5,
        years=[1955, 1956, 1957],
        hindcast_years=[1990, 1991],
        workdir=tmp_path,
    )
    assert hc.sizes == {
        "lead": minimal_emulator.duration + 1,
        "init": 2,
        "member": 3,
        "nbasins": 1,
    }

    # Construct climpred HindcastEnsemble
    qobs = xr.open_dataset(salmon_meteo).qobs
    qobs.name = "q_sim"
    hce = HindcastEnsemble(hc).add_observations(qobs)

    rank_histo_verif = hce.verify(
        metric="rank_histogram",
        comparison="m2o",
        dim=["member", "init"],
        alignment="same_inits",
    )
    assert rank_histo_verif.q_sim.shape[0] == minimal_emulator.duration + 1
