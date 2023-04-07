import datetime as dt

from ravenpy import Emulator, EnsembleReader
from ravenpy.new_config import commands as rc
from ravenpy.new_config import options as o
from ravenpy.new_config.emulators import GR4JCN

# Alternative names for variables in meteo forcing file
alt_names = {
    "RAINFALL": "rain",
    "TEMP_MIN": "tmin",
    "TEMP_MAX": "tmax",
    "PET": "pet",
    "HYDROGRAPH": "qobs",
    "SNOWFALL": "snow",
}


def test_enkf(salmon_meteo, salmon_hru, tmp_path):
    """Test one run of Ensemble Kalman Filter data assimilation."""
    cls = GR4JCN
    # name = "GR4JCN"
    data_type = ["RAINFALL", "TEMP_MIN", "TEMP_MAX", "SNOWFALL"]

    conf = cls(
        params=(0.14, -0.005, 576, 7.0, 1.1, 0.92),
        Gauge=[
            rc.Gauge.from_nc(
                salmon_meteo,
                data_type=data_type,
                alt_names=alt_names,
                data_kwds={"ALL": {"elevation": salmon_hru["land"]["elevation"]}},
            ),
        ],
        ObservationData=[rc.ObservationData.from_nc(salmon_meteo, alt_names="qobs")],
        HRUs=[salmon_hru["land"]],
        StartDate=dt.datetime(1996, 9, 1),
        Duration=30,
        EnsembleMode=rc.EnsembleMode(n=7),
        EnKFMode=o.EnKFMode.SPINUP,
        AssimilationStartTime=dt.datetime(1996, 9, 2),
        RunName="spinup",
        EvaluationMetrics=("NASH_SUTCLIFFE",),
        GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
        OutputDirectoryFormat="./ens_*",
        ForcingPerturbation=[
            rc.ForcingPerturbation(
                forcing="RAINFALL",
                dist="DIST_NORMAL",
                p1=0,
                p2=0.5,
                adj="MULTIPLICATIVE",
            ),
            rc.ForcingPerturbation(
                forcing="SNOWFALL",
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
        ObservationalErrorModel=[
            rc.ObservationalErrorModel(
                state="STREAMFLOW",
                dist="DIST_NORMAL",
                p1=1,
                p2=0.07,
                adj="MULTIPLICATIVE",
            )
        ],
        DebugMode=True,
        NoisyMode=True,
    )

    tmp_path = "/tmp/enkf1"
    conf.zip(tmp_path, overwrite=True)
    # Spin-up

    spinup = Emulator(config=conf, workdir=tmp_path, overwrite=True).run(overwrite=True)

    # Closed Loop
    # loop_conf = conf.set_solution(spinup.solution)
    # conf.enkf_mode = o.EnKFMode.CLOSED_LOOP
    # conf.run_name = "loop"
    # conf.solution_run_name = "spinup"
    # conf.uniform_initial_conditions = None
    # loop = Emulator(config=conf, workdir=tmp_path).run(overwrite=True)

    # Forecast

    # conf.duration =

    # paths = list(tmp_path.glob("ens_*"))
    # paths.sort()
    # assert len(paths) == 7

    # ens = EnsembleReader(conf.run_name, paths)
    # assert len(ens.files["hydrograph"]) == 7
