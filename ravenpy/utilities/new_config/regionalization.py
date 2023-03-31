#!/usr/bin/env python3
"""
Created on Thu Mar 30 12:04:26 2023

@author: ets
"""
import tempfile
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import xarray as xr

from ravenpy.ravenpy import Emulator, EnsembleReader

from .. import coords
from ..regionalization import (
    IDW,
    distance,
    multiple_linear_regression,
    regionalization_params,
    similarity,
)
from .coords import param


def regionalize(
    config: "Config",
    method: str,
    nash,
    params=None,
    props=None,
    target_props=None,
    size: int = 5,
    min_NSE: float = 0.6,
    workdir: Union[str, Path] = None,
    overwrite: bool = False,
    **kwds,
):
    """Perform regionalization for catchment whose outlet is defined by coordinates.

    Parameters
    ----------
    config: Config
      Symbolic emulator configuration. Only GR4JCN, HMETS and Mohyse are supported.
    method : {'MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'}
      Name of the regionalization method to use.
    nash : pd.Series
      NSE values for the parameters of gauged catchments.
    params : pd.DataFrame
      Model parameters of gauged catchments. Needed for all but MRL method.
    props : pd.DataFrame
      Properties of gauged catchments to be analyzed for the regionalization. Needed for MLR and RA methods.
    target_props : pd.Series or dict
      Properties of ungauged catchment. Needed for MLR and RA methods.
    size : int
      Number of catchments to use in the regionalization.
    min_NSE : float
      Minimum calibration NSE value required to be considered as a donor.
    workdir: Union[str, Path]
      Work directory. If None, a temporary directory will be created.
    overwrite: bool
      If True, existing files will be overwritten.
    kwds : {}
      Model configuration parameters, including the forcing files (ts).

    Returns
    -------
    (qsim, ensemble)
    qsim : DataArray (time, )
      Multi-donor averaged predicted streamflow.
    ensemble : Dataset
      q_sim : DataArray  (realization, time)
        Ensemble of members based on number of donors.
      parameter : DataArray (realization, param)
        Parameters used to run the model.
    """
    name = config.__class__.__name__
    if name not in ["GR4JCN", "HMETS", "Mohyse"]:
        raise ValueError(f"Emulator {name} is not supported for regionalization.")

    if not config.is_symbolic:
        raise ValueError("config should be a symbolic configuration.")

    # TODO: Include list of available properties in docstring.
    # TODO: Add error checking for source, target stuff wrt method chosen.
    workdir = Path(workdir or tempfile.mkdtemp())

    # Select properties based on those available in the ungauged properties DataFrame.
    if isinstance(target_props, dict):
        ungauged_properties = pd.Series(target_props)
    elif isinstance(target_props, pd.Series):
        ungauged_properties = target_props
    elif isinstance(target_props, pd.DataFrame):
        ungauged_properties = target_props.to_series()
    else:
        raise ValueError

    cr = coords.realization(1 if method == "MLR" else size)
    cp = param(name)

    # Filter on NSE
    valid = nash > min_NSE
    filtered_params = params.where(valid).dropna()
    filtered_prop = props.where(valid).dropna()

    # Check to see if we have enough data, otherwise raise error
    if len(filtered_prop) < size and method != "MLR":
        raise ValueError(
            "Hydrological_model and minimum NSE threshold \
                         combination is too strict for the number of donor \
                         basins. Please reduce the number of donor basins OR \
                         reduce the minimum NSE threshold."
        )

    # Rank the matrix according to the similarity or distance.
    if method in ["PS", "PS_IDW", "PS_IDW_RA"]:  # Physical similarity
        dist = similarity(filtered_prop, ungauged_properties)
    else:  # Geographical distance.
        dist = distance(filtered_prop, ungauged_properties)

    # Series of distances for the first `size` best donors
    sdist = dist.sort_values().iloc[:size]

    # Pick the donors' model parameters and catchment properties
    sparams = filtered_params.loc[sdist.index]
    sprop = filtered_prop.loc[sdist.index]

    # Get the list of parameters to run
    reg_params = regionalization_params(
        method, sparams, sprop, ungauged_properties, filtered_params, filtered_prop
    )

    # Run the model over all parameters and create ensemble DataArray

    ensemble = []
    for i, rparams in enumerate(reg_params):
        model_config_tmp = config.set_params(rparams)

        out = Emulator(model_config_tmp, workdir=workdir / f"donor_{i}").run(
            overwrite=overwrite
        )

        # Append to the ensemble.
        ensemble.append(out)

    qsim_obj = EnsembleReader(runs=ensemble, dim="members")
    qsims = qsim_obj.hydrograph.q_sim

    # 3. Aggregate runs into a single result -> dataset
    if method in [
        "MLR",
        "SP",
        "PS",
    ]:  # Average (one realization for MLR, so no effect).
        qsim = qsims.mean(dim="members", keep_attrs=True)
    elif (
        "IDW" in method
    ):  # Here we are replacing the mean by the IDW average, keeping attributes and dimensions.
        qsim = IDW(qsims, sdist)
    else:
        raise ValueError(f"No matching algorithm for {method}")

    # Metadata handling
    # TODO: Store the basin_name

    # Create a DataArray for the parameters used in the regionalization
    param_da = xr.DataArray(
        reg_params,
        dims=("members", "param"),
        coords={"param": cp, "members": cr},
        attrs={"long_name": "Model parameters used in the regionalization."},
    )

    ens = xr.Dataset(
        data_vars={"q_sim": qsims, "parameter": param_da},
        attrs={
            "title": "Regionalization ensemble",
            # "institution": "",
            "source": f"Hydrological model {name}",
            "history": "Created by ravenpy regionalize.",
            # "references": "",
            "comment": f"Regionalization method: {method}",
        },
    )

    # TODO: Add global attributes (model name, date, version, etc)
    return qsim, ens
