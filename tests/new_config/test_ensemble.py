import datetime as dt

from ravenpy import Emulator, EnsembleReader
from ravenpy.new_config import commands as rc
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
        Gauge=rc.Gauge.from_nc(
            salmon_meteo,
            data_type=data_type,
            alt_names=alt_names,
            extra={1: {"elevation": salmon_hru["land"]["elevation"]}},
        ),
        ObservationData=rc.ObservationData.from_nc(salmon_meteo, alt_names="qobs"),
        HRUs=[salmon_hru["land"]],
        StartDate=dt.datetime(1996, 9, 1),
        EndDate=dt.datetime(1996, 9, 15),
        AssimilationStartTime=dt.datetime(1996, 9, 8),
        RunName="test",
        EvaluationMetrics=("NASH_SUTCLIFFE",),
        GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
        EnsembleMode=rc.EnsembleMode(n=7),
        OutputDirectoryFormat=tmp_path / "ens_*",
        ForcingPerturbation=[
            rc.ForcingPerturbation(
                forcing="RAINFALL",
                dist="DIST_NORMAL",
                p1=1,
                p2=0.05,
                adj="MULTIPLICATIVE",
            ),
            rc.ForcingPerturbation(
                forcing="SNOWFALL",
                dist="DIST_NORMAL",
                p1=1,
                p2=0.05,
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

    e = Emulator(config=conf, workdir=tmp_path)
    e.write_rv()
    e.run()
    paths = list(tmp_path.glob("ens_*"))
    paths.sort()
    assert len(paths) == 14

    ens = EnsembleReader(conf.run_name, paths)
    assert len(ens.files["hydrograph"]) == 14
