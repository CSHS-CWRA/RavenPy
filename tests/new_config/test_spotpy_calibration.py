import datetime as dt
import time

import spotpy

from ravenpy.new_config import commands as rc
from ravenpy.new_config.emulators.gr4jcn import GR4JCN
from ravenpy.new_config.emulators.hmets import HMETS
from ravenpy.new_config.emulators.mohyse import Mohyse
from ravenpy.utilities.new_config.calibration import SpotSetup

salmon_river = "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"


class TestGR4JCNSpotpy:
    def test_simple(self, get_file, tmpdir):
        ts = get_file(salmon_river)

        alt_names = {
            "RAINFALL": "rain",
            "TEMP_MIN": "tmin",
            "TEMP_MAX": "tmax",
            "PET": "pet",
            "HYDROGRAPH": "qobs",
            "SNOWFALL": "snow",
        }

        salmon_land_hru_1 = dict(
            area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659
        )

        model = GR4JCN(
            # params=[0.529, -3.396, 407.29, 1.072, 16.9, 0.947],  # Pas d'erreur?
            ObservationData=[rc.ObservationData.from_nc(
                ts, alt_names="qobs"
            )],  # Juste necessaire pour la calibration
            Gauge=[rc.Gauge.from_nc(
                ts,
                alt_names=alt_names,
                extra={1: {"elevation": salmon_land_hru_1["elevation"]}},
            )],
            HRUs=[salmon_land_hru_1],
            StartDate=dt.datetime(1960, 1, 1),
            EndDate=dt.datetime(2002, 1, 1),
            RunName="test",  # Le nom par defaut?
            EvaluationMetrics=("NASH_SUTCLIFFE",),  # Devrait pas etre obligatoire
        )

        spot_setup = SpotSetup(
            config=model,
            low=(0.01, -15.0, 10.0, 0.0, 1.0, 0.0),
            high=(2.5, 10.0, 700.0, 7.0, 30.0, 1.0),
        )

        sampler = spotpy.algorithms.dds(
            spot_setup, dbname="RAVEN_model_run", dbformat="ram", save_sim=False
        )
        rep = 8

        sampler.sample(rep, trials=1)        
        assert spot_setup.diagnostics is not None
        
        results = sampler.getdata()
        assert len(spotpy.analyser.get_best_parameterset(results)[0])==6
