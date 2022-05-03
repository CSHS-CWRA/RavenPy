import datetime as dt
import time

import spotpy

from ravenpy.models import GR4JCN
from ravenpy.utilities.calibration import SpotpySetup
from ravenpy.utilities.testdata import get_local_testdata


class TestGR4JCN_Spotpy:
    def test_simple(self):

        model = GR4JCN()

        TS = get_local_testdata(
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
        )

        salmon_land_hru_1 = dict(
            area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659
        )

        model.config.rvh.hrus = (GR4JCN.LandHRU(**salmon_land_hru_1),)

        # Parameter bounds
        model.low = (0.01, -15.0, 10.0, 0.0, 1.0, 0.0)
        model.high = (2.5, 10.0, 700.0, 7.0, 30.0, 1.0)

        model.config.rvi.start_date = dt.datetime(2000, 1, 1)
        model.config.rvi.end_date = dt.datetime(2002, 1, 1)
        model.config.rvi.run_name = "test"

        spot_setup = SpotpySetup(model=model, ts=TS, obj_func=None)
        sampler = spotpy.algorithms.dds(
            spot_setup, dbname="RAVEN_model_run", dbformat="ram", save_sim=False
        )
        rep = 100

        tic = time.time()
        sampler.sample(rep, trials=1)
        toc = time.time()
        total_time = toc - tic
        print(total_time)
        # 10 evals = 4.61 seconds
        # 100 evals = 49 seconds
        # 1000 evals = 763 seconds
        # 1200 evals = 1000 seconds
        # 2500 evals = 2951 seconds
        # 5000 evals = 9580 seconds
        # 7500 evals = 20004 seconds
        # 10000 evals = 34990 seconds

        results = sampler.getdata()
        model.config.update("params", spotpy.analyser.get_best_parameterset(results)[0])
        model(TS)
        objfun = model.diagnostics["DIAG_NASH_SUTCLIFFE"][0]
        print(objfun)
