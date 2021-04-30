import datetime as dt

from ravenpy.models import GR4JCN
from ravenpy.utilities.testdata import get_local_testdata

"""
Test to perform a hindcast using auto-queried ECCC data aggregated on THREDDS.
Currently only runs GEPS, eventually will run GEPS, GDPS, REPS and RDPS.
To do so will need to add the actual data from ECCC but this is a proof of concept.
"""

hru = GR4JCN.LandHRU(
    area=44250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659
)


class TestECCCForecast:
    def test_forecasting_GEPS(self):

        # Prepare a RAVEN model run using historical data, GR4JCN in this case.
        # This is a dummy run to get initial states. In a real forecast situation,
        # this run would end on the day before the forecast, but process is the same.
        ts = get_local_testdata(
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
        )
        model = GR4JCN()
        model(
            ts,
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2002, 6, 1),
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
            hrus=(hru,),
        )

        # Extract the final states that will be used as the next initial states
        rvc = model.outputs["solution"]

        # Collect test forecast data for location and climate model (20 members)
        ts20 = get_local_testdata("eccc_forecasts/geps_watershed.nc")
        nm = 20

        # It is necessary to clean the model state because the input variables of the previous
        # model are not the same as the ones provided in the forecast model. therefore, if we
        # do not clean, the model will simply add the hindcast file to the list of available
        # data provided in the testdata above. Then the dates will not work, and the model errors.

        model = GR4JCN()

        model.config.rvc.parse_solution(rvc.read_text())

        model(
            ts=(ts20,),
            nc_index=range(nm),
            duration=9,
            params=(0.529, -3.396, 407.29, 1.072, 16.9, 0.947),
            hrus=(hru,),
            overwrite=True,
            pr={"time_shift": -0.25, "deaccumulate": True},
            tas={"time_shift": -0.25},
        )

        # The model now has the forecast data generated and it has 10 days of forecasts.
        assert len(model.q_sim.values) == 10

        # Also see if GEPS has 20 members produced.
        assert model.q_sim.values.shape[1] == nm

        # Check all members are different (checking snow because data in winter)
        assert len(set(model.storage.Snow.isel(time=-1).values)) == nm
