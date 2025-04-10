"""Tools for hydrological regionalization."""

import logging
import tempfile
from pathlib import Path
from typing import Any, Callable, Optional, Union

import numpy as np
import pandas as pd
import statsmodels.api as sm
import xarray as xr
from haversine import haversine_vector

from ravenpy.config import Config
from ravenpy.ravenpy import Emulator, EnsembleReader

from . import coords

# from statsmodels.regression.linear_model import RegressionResults


LOGGER = logging.getLogger("PYWPS")

regionalisation_data_dir = Path(__file__).parent.parent / "data" / "regionalisation"


def regionalize(
    config: Config,
    method: str,
    nash: pd.Series,
    params: pd.DataFrame = None,
    props: pd.DataFrame = None,
    target_props: Union[pd.Series, dict] = None,
    size: int = 5,
    min_NSE: float = 0.6,  # noqa: N803
    workdir: Optional[Union[str, Path]] = None,
    overwrite: bool = False,
) -> tuple[xr.DataArray, xr.Dataset]:
    r"""Perform regionalization for catchment whose outlet is defined by coordinates.

    Parameters
    ----------
    config : ravenpy.config.rvs.Config
        Symbolic emulator configuration. Only GR4JCN, HMETS and Mohyse are supported.
    method : {'MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'}
        Name of the regionalization method to use.
    nash : pd.Series
        NSE values for the parameters of gauged catchments.
    params : pd.DataFrame
        Model parameters of gauged catchments. Needed for all but MRL method.
    props : pd.DataFrame
        Properties of gauged catchments to be analyzed for the regionalization.
        Needed for MLR and RA methods.
    target_props : pd.Series or dict
        Properties of ungauged catchment. Needed for MLR and RA methods.
    size : int
        Number of catchments to use in the regionalization.
    min_NSE : float
        Minimum calibration NSE value required to be considered as a donor.
    workdir : Union[str, Path]
        Work directory. If None, a temporary directory will be created.
    overwrite : bool
        If True, existing files will be overwritten.

    Returns
    -------
    qsim : DataArray (time, )
        Multi-donor averaged predicted streamflow.
    ensemble : Dataset
        A Dataset containing the ensemble of simulations and parameters used:

        - q_sim : DataArray (realization, time)
          Ensemble of members based on number of donors.

        - parameter : DataArray (realization, param)
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
    cp = coords.param(name)

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


def read_gauged_properties(properties) -> pd.DataFrame:
    """Return table of gauged catchments properties over North America.

    Returns
    -------
    pd.DataFrame
        Catchment properties keyed by catchment ID.
    """
    f = regionalisation_data_dir / "gauged_catchment_properties.csv"
    proptable = pd.read_csv(f, index_col="ID")

    return proptable[properties]


def read_gauged_params(model):
    """Return table of NASH-Sutcliffe Efficiency values and model parameters for North American catchments.

    Returns
    -------
    pd.DataFrame
        Nash-Sutcliffe Efficiency keyed by catchment ID.
    pd.DataFrame
        Model parameters keyed by catchment ID.
    """
    f = regionalisation_data_dir / f"{model}_parameters.csv"
    params = pd.read_csv(f, index_col="ID")

    return params["NASH"], params.iloc[:, 1:]


def distance(gauged: pd.DataFrame, ungauged: pd.Series) -> pd.Series:
    """Return geographic distance [km] between ungauged and database of gauged catchments.

    Parameters
    ----------
    gauged : pd.DataFrame
        Table containing columns for longitude and latitude of catchment's centroid.
    ungauged : pd.Series
        Coordinates of the ungauged catchment.

    Returns
    -------
    pd.Series
    """
    gauged_array = np.array(list(zip(gauged.latitude.values, gauged.longitude.values)))

    return pd.Series(
        data=haversine_vector(
            gauged_array, np.array([ungauged.latitude, ungauged.longitude]), comb=True
        )[0],
        index=gauged.index,
    )


def similarity(
    gauged: pd.DataFrame, ungauged: pd.DataFrame, kind: str = "ptp"
) -> pd.Series:
    """Return similarity measure between gauged and ungauged catchments.

    Parameters
    ----------
    gauged : pd.DataFrame
        Gauged catchment properties.
    ungauged : pd.DataFrame
        Ungauged catchment properties.
    kind : {'ptp', 'std', 'iqr'}
        Normalization method: peak to peak (maximum - minimum), standard deviation, inter-quartile range.

    Returns
    -------
    pd.Series
    """
    stats = gauged.describe()

    if kind == "ptp":
        spread = stats.loc["max"] - stats.loc["min"]
    elif kind == "std":
        spread = stats.loc["std"]
    elif kind == "iqr":
        spread = stats.loc["75%"] - stats.loc["25%"]

    d = ungauged.values - gauged.values
    n = np.abs(d) / spread.values
    return pd.Series(data=n.sum(axis=1), index=gauged.index)


# FIXME: gauged_properties is not used
def regionalization_params(
    method: str,
    gauged_params: pd.DataFrame,
    gauged_properties: pd.DataFrame,  # noqa: F841
    ungauged_properties: pd.DataFrame,
    filtered_params: pd.DataFrame,
    filtered_prop: pd.DataFrame,
) -> Union[list[list[float]], list[np.ndarray]]:
    """
    Return the model parameters to use for the regionalization.

    Parameters
    ----------
    method : {'MLR', 'SP', 'PS', 'SP_IDW', 'PS_IDW', 'SP_IDW_RA', 'PS_IDW_RA'}
        Name of the regionalization method to use.
    gauged_params : pd.DataFrame
        A DataFrame of parameters for donor catchments (size = number of donors).
    gauged_properties : pd.DataFrame
        A DataFrame of properties of the donor catchments  (size = number of donors).
    ungauged_properties : pd.DataFrame
        A DataFrame of properties of the ungauged catchment (size = 1).
    filtered_params : pd.DataFrame
        A DataFrame of parameters of all filtered catchments (size = all catchments with NSE > min_NSE).
    filtered_prop : pd.DataFrame
        A DataFrame of properties of all filtered catchments (size = all catchments with NSE > min_NSE).

    Returns
    -------
    list
        A list of model parameters to be used for the regionalization.
    """
    if method == "MLR" or "RA" in method:
        mlr_params, r2 = multiple_linear_regression(
            filtered_prop, filtered_params, ungauged_properties.to_frame().T
        )

        if method == "MLR":  # Return the multiple linear regression parameters.
            out = [
                mlr_params,
            ]

        else:
            gp = gauged_params.copy()

            for p, r, col in zip(mlr_params, r2, gauged_params):
                # If we have an R2 > 0.5 then we consider this to be a better estimator

                if r > 0.5:
                    gp[col] = p

            out = gp.values

    else:
        out = gauged_params.values

    return out


def IDW(qsims: xr.DataArray, dist: pd.Series) -> xr.DataArray:  # noqa: N802
    """Inverse distance weighting.

    Parameters
    ----------
    qsims : xr.DataArray
        Ensemble of hydrogram stacked along the `members` dimension.
    dist : pd.Series
        Distance from catchment which generated each hydrogram to target catchment.

    Returns
    -------
    xr.DataArray
        Inverse distance weighted average of ensemble.
    """
    # In IDW, weights are 1 / distance
    weights = xr.DataArray(
        1.0 / dist, dims="members", coords={"members": qsims.members}
    )

    # Make weights sum to one
    weights /= weights.sum(axis=0)

    # Calculate weighted average.
    out = qsims.dot(weights)
    out.name = qsims.name
    out.attrs = qsims.attrs
    return out


def multiple_linear_regression(
    source: pd.DataFrame, params: pd.DataFrame, target: pd.DataFrame
) -> tuple[list[Any], list[Callable[[], Any]]]:
    """Multiple Linear Regression for model parameters over catchment properties.

    Uses known catchment properties and model parameters to estimate model parameter over an
    ungauged catchment using its properties.

    Parameters
    ----------
    source : pd.DataFrame
        Properties of gauged catchments.
    params : pd.DataFrame
        Model parameters of gauged catchments.
    target : pd.DataFrame
        Properties of the ungauged catchment.

    Returns
    -------
    list of Any, list of Callable or Any
        A named tuple of the estimated model parameters and the R2 of the linear regression.
    """
    # Add constants to the gauged predictors
    x = sm.add_constant(source)

    # Add the constant '1' for the ungauged catchment predictors
    predictors = sm.add_constant(target, prepend=True, has_constant="add")

    # Perform regression for each parameter
    regression = [sm.OLS(params[param].values, x).fit() for param in params]

    # Perform prediction on each parameter based on the predictors
    mlr_parameters = [r.predict(exog=predictors)[0] for r in regression]

    # Extract the adjusted r_squared value for each parameter
    r2 = [r.rsquared_adj for r in regression]

    return mlr_parameters, r2
