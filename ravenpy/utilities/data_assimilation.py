"""
Created on Sun Nov  1 20:48:03 2020

@author: Richard Arsenault
"""
import datetime as dt
import math
import os
from copy import deepcopy
from dataclasses import replace
from typing import Dict, List, Sequence, Tuple, Union

import numpy as np
import xarray as xr

"""
model = Raven model instance, preset with parameters etc.
xa = the set of state variables. In this case, soil0 and soil1 from GR4JCN. Will eventually need to update for other models, other variables
ts = the timeseries of inputs needed by model (tas, pr, qobs, etc.)
days = data for the days of assimilation. Assimilation is performed after last day.
number_members = number of EnKF members
std = variables for uncertainty estimation of input and streamflow variables.
    ex: std = {"rainfall": 0.30,"prsn": 0.30, "tasmin": 2.0, "tasmax": 2.0, "water_volume_transport_in_river_channel": 0.15}
precip = standard deviation used to sample precip, uses gamma distribution (fraction of observed value)
temperature = standard deviation used to sample temperature, normal dist (degrees Celcius)
qobs = standard deviation used to sample observed streamflow, normal dist (fraction of observed value)
"""


def assimilate(
    model,
    ts: Union[str, os.PathLike],
    q_obs: xr.Dataset,
    keys: Tuple,
    basin_states: Sequence,
    hru_states: Sequence,
    days: dt.datetime,
):
    """Assimilate streamflow over one day.

    Parameters
    ----------
    model : raven.Model
      Raven model instance configured to run.
    ts : str, Path
      Perturbed time series.
    keys : tuple
      Name of hru_state attributes to be assimilated, for example ("soil0", "soil1").
    basin_states : sequence
      Model initial conditions, BasinStateVariables instances.
    hru_states : sequence
      Model initial conditions, HRUStateVariables instances.
    q_obs : xarray.Dataset
      The actual observed streamflow over the entire period.
    ts : str, Path, list
      Input netCDF file names.
    days : datetime.datetime
      Dates ...
    """
    qkey = "water_volume_transport_in_river_channel"

    if len(basin_states) != len(hru_states):
        raise ValueError("`basin_states` and `hru_states` must have the same length.")

    model = deepcopy(model)

    # Number of members
    n_members = len(basin_states)

    # Run simulation with perturbed inputs
    model(
        ts,
        parallel=dict(
            hru_state=hru_states, basin_state=basin_states, nc_index=range(n_members)
        ),
    )

    # Extract final states (n_states, n_members)
    f_hru_states, f_basin_states = model.get_final_state()
    x_matrix = np.array(
        [[getattr(state, key) for key in keys] for state in f_hru_states]
    ).T

    # Sanity check
    if x_matrix.shape != (len(keys), n_members):
        raise ValueError

    # vnames = [HRU_NC_MAP[k] for k in keys]
    # x_matrix = model.storage.isel(time=-1)[vnames].to_array()

    # Last time step for assimilation
    # The parallel dimension is state, not nbasins.
    qsim_vector = model.q_sim.isel(time=-1, nbasins=0).values

    # If there are problems related to missing Qobs or other variables, do not assimilate.
    with xr.open_dataset(ts) as perturbed:
        if perturbed.isnull().any().to_array().any():
            xa = x_matrix
        else:
            # Prepare the perturbed Qobs and Qobs errors for each member
            qobs_pert = perturbed[qkey].sel(time=days[-1])
            qobs_error = q_obs.sel(time=days[-1]) - qobs_pert

            # Do the assimilation, return the assimilated states for the next period
            xa = update_state(
                x_matrix, qobs_pert.values, qobs_error.values, qsim_vector
            )

    return [xa, model]


def perturbation(
    da: xr.DataArray, dist: str, std: float, seed: int = None, **kwargs: Dict
):
    """Return perturbed time series.

    Parameters
    ----------
    da : xr.DataArray
      Input time series to be perturbed.
    dist : {"norm", "gamma", "rnorm"}
      Name of statistical distribution from which random perturbations are drawn. `rnorm` stands for relative normal, where the given standard deviation is multiplied by the value being perturbed.
    std : float
      Standard deviation of random perturbation.
    seed : int
      Seed for the random number generator. Setting the same seed for different variables will ensure the same
      perturbations are generated.
    kwargs : dict
      Name and size of additional dimensions, apart from time.

    """
    nt = len(da.time)
    size = list(kwargs.values()) + [nt]
    dims = list(kwargs.keys()) + ["time"]

    # Create random state object. If seed is None, the seed itself will be random.
    rs = np.random.RandomState(seed)

    if dist == "norm":
        r = rs.normal(0, std, size=size)
        out = da + xr.DataArray(r, dims=dims, coords={"time": da.time})

    elif dist == "rnorm":
        r = rs.normal(0, da * std, size=size)
        out = da + xr.DataArray(r, dims=dims, coords={"time": da.time})

    elif dist == "gamma":
        shape = (da**2) / (std * da) ** 2
        scale = ((std * da) ** 2) / da
        r = np.nan_to_num(rs.gamma(shape=shape, scale=scale, size=size), nan=0.0)
        out = xr.DataArray(r, dims=dims, coords={"time": da.time})

    else:
        raise AttributeError(f"{dist} is not supported.")

    out.attrs.update(da.attrs)
    return out


def update_state(
    x: np.ndarray, qobs_pert: np.ndarray, qobs_error: np.ndarray, qsim: np.ndarray
):
    """
    Update model state by assimilation.

    Parameters
    ----------
    x : ndarray (n_states, n_members)
      Model state initial values.
    qobs_pert : ndarray (n_members)
      Perturbed observed streamflows.
    qobs_error : ndarray (n_members)
      Perturbation added to qobs to get qobs_pert.
    qsim : ndarray (n_members)
      Simulated streamflows.

    Returns
    -------
    ndarray (n_states, n_members)
      Model state values after assimilation.

    References
    ----------
    The Ensemble Kalman Filter: theoretical formulation and practical implementation, Evensen 2003
    https://link.springer.com/article/10.1007%2Fs10236-003-0036-9

    University of colorado report on the efficient implementation of EnKF, Jan Mandel, 2006
    http://ccm.ucdenver.edu/reports/rep231.pdf
    """
    n_states, n_members = np.shape(x)

    # Make sure arrays have shape (1, n_members)
    qobs_pert = np.atleast_2d(qobs_pert)
    qobs_error = np.atleast_2d(qobs_error)
    qsim = np.atleast_2d(qsim)

    z = np.dot(qsim, np.ones((n_members, 1)))
    ha = qsim - (z * np.ones((1, n_members))) / n_members
    y = qobs_pert - qsim

    # Equations 4.1 from Mandel, 2006
    re = np.dot(qobs_error, qobs_error.transpose()) / n_members
    p = re + (np.dot(ha, ha.transpose())) / (n_members - 1)
    m = np.dot(p**-1, y)
    z = (ha.transpose()) * m
    a = (
        x
        - np.dot((np.dot(x, np.ones((n_members, 1)))), np.ones((1, n_members)))
        / n_members
    )
    xa = x + (np.dot(a, z)) / (n_members - 1)
    xa = np.maximum(xa, 0)

    return xa


def perturb_full_series(
    model,
    std: Dict,
    start_date: dt.datetime,
    end_date: dt.datetime,
    dists: Tuple,
    n_members: int = 25,
):
    """
    Create an ensemble of randomly perturbed input timeseries for the given model.
    Each input variable is modified by 'n_members' time series of random values
    drawn from a given distribution and parameters.

    Parameters
    ----------
    model : ravenpy.Raven instance
      The model that will be used to perform the simulations and assimilation
    std : dict
      Standard deviation of the perturbation noise for each input variable, keyed by variable standard name.
    start_date : datetime.datetime
      Start date of the assimilation run.
    end_date : datetime.datetime
      End date of the assimilation run.
    dists : tuple
      Pairs of variable:distribution to identify the type of distribution each variable follows.
    n_members : int
      Number of members to use in the Ensemble Kalman Filter.

    Returns
    -------
    xarray.Dataset
      Dataset containing all perturbed hydrometeorological timeseries for the entire period.
    """

    # Use the same random seed for both tasmin and tasmax
    rs = np.random.SeedSequence(None).generate_state(1)[0]
    seed = {"tasmin": rs, "tasmax": rs}

    # === Create perturbed time series for full assimilation period ====
    perturbed = {}
    for key, s in std.items():
        nc = model.config.rvt._var_cmds[key]

        with xr.open_dataset(nc.file_name_nc) as ds:
            da = ds.get(nc.var_name_nc).sel(time=slice(start_date, end_date))

            perturbed[key] = perturbation(
                da,
                dists.get(key, "norm"),
                std=s,
                seed=seed.get(key, None),
                member=n_members,
            )

    return perturbed


def assimilation_initialization(
    model,
    ts: str,
    start_date: dt.datetime,
    end_date: dt.datetime,
    area: float,
    elevation: float,
    latitude: float,
    longitude: float,
    params: List[float],
    assim_var: List[str],
    n_members: int = 25,
):
    """
    This function takes a Raven model setup and runs it for a given period of time
    as defined by start_date and end_date. At the end of the period, the initial
    states are duplicated n_members times and returned so that the next
    assimilation step can have access to an ensemble of initial states.

    Parameters
    ----------
    model : ravenpy.Raven instance
      The model that will be used to perform the simulations and assimilation.
    ts : str
      Path to the forcing data timeseries. For assimilation, this is perturbed data.
    start_date : datetime.datetime
      Start date of the period used to initialize states for assimilation.
    end_date : datetime.datetime
      Date at which we will generate states for the upcoming assimilation.
    area : float
      Catchment area, in square kilometers (km^2)
    elevation : float
      Catchment elevation, in meters (m).
    latitude : float
      Catchment latitude, in degrees.
    longitude : float
      Catchment longitude, in degrees (negative west).
    params : List[float]
      The hydrological model parameters used for the assimilation and simulation.
    assim_var : List[str]
      List of hydrological model internal variables to be modified during assimilation.
    n_members : int
      Number of ensemble members for the Ensemble Kalman Filter. The default value is 25.

    Returns
    -------
    ravenpy.Raven instance
      The hydrological model with the internal states after the n_members simulations.
    np.array
      Array of state variables used overwrite the current initial states.
    ravenpy.Raven.hru_states instance
      The Raven model states for the hru information (size n_members)
    ravenpy.Raven.basin_states instance
      The Raven model states for the basin information (size n_members).
    """
    # Set model options
    model(
        ts=ts,
        start_date=start_date,
        end_date=end_date,
        area=area,
        elevation=elevation,
        latitude=latitude,
        longitude=longitude,
        params=params,
    )
    """
    # This section of code was the old way of doing it. Should be removed
    # when notebooks will be working with the new implementations.

    model.config.rvh.hrus = (
        GR4JCN.LandHRU(
            area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659
        ),
    )

    model.config.rvp.params = GR4JCN.Params(
        0.1353389, -0.005067198, 576.8007, 6.986121, 1.102917, 0.9224778
    )  # SALMON

    # ==== Initialization (just to get reasonable states) ====
    # Set initialization run options
    model.config.rvi.run_name = "init"
    model.config.rvi.start_date = start_date
    # Richard: what is the end date policy for the init run ?
    model.config.rvi.end_date = start_date + dt.timedelta(days=assim_days[0])

    # Run the model
    model([ts])
    """
    # Extract final model states
    hru_state, basin_state = model.get_final_state()
    xa = n_members * [getattr(hru_state, key) for key in assim_var]
    hru_states = n_members * [hru_state]
    basin_states = n_members * [basin_state]

    return model, xa, hru_states, basin_states


def sequential_assimilation(
    model,
    hru_states,
    basin_states,
    p_fn,
    q_obs,
    assim_var,
    start_date,
    end_date,
    n_members=25,
    assim_step_days=7,
):
    """
    This function performs data assimilation at regular intervals over a long
    period, as defined by start_date and end_date. Intervals are equal to
    assim_step_days and assimilation is performed every assim_step_days until
    the end of the period. The function returns the complete assimilated
    streamflow for each of the n_members desired in the Ensemble Kalman Filter.

    Parameters
    ----------
    model : ravenpy.Raven instance
      The hydrological model with the internal states after the n_members simulations.
    hru_states : ravenpy.Raven.hru_states instance
      The Raven model states containing the initial ensemble hru information (size n_members)
    basin_states : ravenpy.Raven.basin_states instance
      The Raven model states containing the initial ensemble basin information (size n_members).
    p_fn : string
      Path to the perturbed forcing data netcdf file.
    q_obs : xarray.Dataset
      The actual observed streamflow over the entire period.
    assim_var : tuple of strings
      List of hydrological model internal variables to be modified during assimilation.
    start_date : datetime.datetime
      Start date of the sequence over which we want to perform regular data assimilation.
    end_date : datetime.datetime
      Date of the end of the sequence over which we want to perform regular data assimilation.
    n_members : integer, optional
      Number of ensemble members for the Ensemble Kalman Filter. The default is 25.
    assim_step_days : integer, optional
      Number of days between each data assimilation step. Days in between are simply simulated without assimilation.
      The default value is 7.

    Returns
    -------
    xarray.DataArray
        Array of assimilated streamflows for the full period duration. Size is n_members x time.
    ravenpy.Raven.hru_states instance
        The Raven model states for the hru information at the end of the period (size n_members)
    ravenpy.Raven.basin_states instance
        The Raven model states for the basin information at the end of the period (size n_members).
    """

    # ==== Assimilation ====
    q_assim = []
    sd = start_date

    # Make list of days on which we will assimilate
    number_days = end_date - start_date
    number_days = number_days.days
    assim_dates = math.floor(number_days / assim_step_days) * [assim_step_days]

    for i, ndays in enumerate(assim_dates):

        dates = [sd + dt.timedelta(days=x) for x in range(ndays)]
        model.config.rvi.end_date = dates[-1]
        model.config.rvi.run_name = f"assim_{i}"

        # Perform the first assimilation step here
        [xa, model] = assimilate(
            model, p_fn, q_obs, assim_var, basin_states, hru_states, dates
        )

        # Save streamflow simulation
        q_assim.append(model.q_sim.isel(nbasins=0))

        # Update the start-time for the next loop
        sd += dt.timedelta(days=ndays)
        model.config.rvi.start_date = sd

        # Get new initial conditions and feed assimilated values
        hru_states, basin_states = model.get_final_state()
        hru_states = [
            replace(hru_states[i], **dict(zip(assim_var, xa[:, i])))
            for i in range(n_members)
        ]

    q_assim = xr.concat(q_assim, dim="time")

    return q_assim, hru_states, basin_states
