import datetime as dt

from ravenpy import Emulator, EnsembleReader
from ravenpy.config import commands as rc
from ravenpy.config import options as o
from ravenpy.config.emulators import GR4JCN

# Alternative names for variables in meteo forcing file
alt_names = {
    "RAINFALL": "rain",
    "TEMP_MIN": "tmin",
    "TEMP_MAX": "tmax",
    "SNOWFALL": "snow",
}


# @pytest.mark.xfail()
def test_enkf(salmon_hru, tmp_path, yangtze):
    """Test one run of Ensemble Kalman Filter data assimilation."""
    cls = GR4JCN
    # name = "GR4JCN"
    data_type = ["RAINFALL", "TEMP_MIN", "TEMP_MAX", "SNOWFALL"]

    salmon_file = yangtze.fetch(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )

    conf = cls(
        params=(0.14, -0.005, 576, 7.0, 1.1, 0.92),
        Gauge=[
            rc.Gauge.from_nc(
                salmon_file,
                data_type=data_type,
                alt_names=alt_names,
                data_kwds={"ALL": {"elevation": salmon_hru["land"]["elevation"]}},
            ),
        ],
        ObservationData=[rc.ObservationData.from_nc(salmon_file, alt_names="qobs")],
        HRUs=[salmon_hru["land"]],
        StartDate=dt.datetime(1996, 9, 1),
        EndDate=dt.datetime(1996, 9, 30),
        EnsembleMode=rc.EnsembleMode(n=7),
        EnKFMode=o.EnKFMode.SPINUP,
        RunName="spinup",
        EvaluationMetrics=("NASH_SUTCLIFFE",),
        OutputDirectoryFormat="./ens_*",
        ForcingPerturbation=[
            rc.ForcingPerturbation(
                forcing="PRECIP",
                dist="DIST_NORMAL",
                p1=1,
                p2=0.07,
                adj="MULTIPLICATIVE",
            ),
        ],
        DefineHRUGroups=["All"],
        HRUGroup=[{"name": "All", "groups": ["1"]}],
        AssimilatedState=[
            rc.AssimilatedState(state="SOIL[0]", group="All"),
            rc.AssimilatedState(state="SOIL[1]", group="All"),
        ],
        AssimilateStreamflow=[rc.AssimilateStreamflow(sb_id=1)],
        ObservationErrorModel=[
            rc.ObservationErrorModel(
                state="STREAMFLOW",
                dist="DIST_NORMAL",
                p1=1,
                p2=0.07,
                adj="MULTIPLICATIVE",
            )
        ],
    )

    # Spinup
    Emulator(config=conf, workdir=tmp_path, overwrite=True).run(overwrite=True)

    # Closed Loop
    conf_loop = conf.duplicate(
        EnKFMode=o.EnKFMode.CLOSED_LOOP,
        RunName="loop",
        SolutionRunName="spinup",
        UniformInitialConditions=None,
        StartDate=dt.datetime(1996, 9, 30),
        EndDate=dt.datetime(1996, 10, 15),
    )

    # Loop
    Emulator(config=conf_loop, workdir=tmp_path).run()

    # Forecast
    conf_cast = conf_loop.duplicate(
        EnKFMode=o.EnKFMode.FORECAST, RunName="forecast", SolutionRunName="loop"
    )

    # Cast
    Emulator(config=conf_cast, workdir=tmp_path).run()

    paths = list(tmp_path.glob("ens_*"))
    paths.sort()
    assert len(paths) == 7

    ens = EnsembleReader(run_name=conf_cast.run_name, paths=paths)
    assert conf_cast.run_name == "forecast"
    assert len(ens.files["hydrograph"]) == 7
