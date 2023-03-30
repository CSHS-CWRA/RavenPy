import datetime as dt

from ravenpy import Emulator, EnsembleReader
from ravenpy.new_config import commands as rc
from ravenpy.new_config.emulators import GR4JCN
from ravenpy.utilities.new_config.forecasting import warm_up

"""
Test to perform a hindcast using Caspar data on THREDDS.
Currently only runs GEPS, eventually will run GEPS, GDPS, REPS and RDPS.
To do so will need to add the actual data from Caspar, currently being downloaded
but this is a good proof of concept.
"""


class TestHindcasting:
    def test_hindcasting_GEPS(self, get_local_testdata, salmon_hru, tmp_path):
        ts20 = get_local_testdata("caspar_eccc_hindcasts/geps_watershed.nc")

        data_type = [
            "RAINFALL",
            "TEMP_MIN",
            "TEMP_MAX",
            "SNOWFALL",
        ]
        alt_names = {
            "RAINFALL": "rain",
            "TEMP_MIN": "tmin",
            "TEMP_MAX": "tmax",
            "PET": "pet",
            "HYDROGRAPH": "qobs",
            "SNOWFALL": "snow",
        }
        hru = salmon_hru["land"]
        data_kwds = {
            "ALL": {"Latitude": hru["latitude"], "Longitude": hru["longitude"]},
            "PRECIP": {
                "Deaccumulate": True,
                "TimeShift": -0.25,
                "LinearTransform": {"scale": 1000},
            },
        }

        pre_conf = GR4JCN(
            StartDate=dt.datetime(2018, 6, 10),
            EndDate=dt.datetime(2018, 6, 14),
            HRUs=[
                hru,
            ],
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
            GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
            Gauge=[
                rc.Gauge.from_nc(
                    ts20, data_type=["PRECIP", "TEMP_AVE"], data_kwds=data_kwds
                ),
            ],
        )

        init_conf = warm_up(pre_conf, duration=8, path=tmp_path / "wup")

        # Create gauges for each member and run the model
        nm = 3
        ens = []
        for i in range(nm):
            g = rc.Gauge.from_nc(
                ts20,
                station_idx=i + 1,
                data_type=["PRECIP", "TEMP_AVE"],
                data_kwds=data_kwds,
            )
            assert g.data[0].read_from_netcdf.station_idx == i + 1
            assert g.data[0].read_from_netcdf.linear_transform.scale == 1000

            mem_conf = init_conf.copy(
                update=dict(
                    gauge=[g],
                    custom_output={
                        "time_per": "DAILY",
                        "stat": "AVERAGE",
                        "variable": "PRECIP",
                        "space_agg": "BY_BASIN",
                    },
                ),
            )

            e = Emulator(
                config=mem_conf,
                workdir=tmp_path / f"m{i:02}",
            )
            ens.append(e.run())

        out = EnsembleReader(runs=ens)

        # The model now has the forecast data generated and it has 5 days of forecasts.
        assert len(out.hydrograph.member) == 3
        assert len(out.hydrograph.time) == 5
