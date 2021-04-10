#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 22:25:32 2021

@author: ets
"""

import datetime as dt
import logging

from ravenpy.models.emulators import GR4JCN
from ravenpy.utilities.forecasting import (
    make_climpred_hindcast_object,
    make_ESP_hindcast_dataset,
)
from ravenpy.utilities.testdata import get_local_testdata


class TestClimpredHindcastVerification:
    def test_simple(self):

        # We don't want climpred logging in the test output
        logger = logging.getLogger()
        logger.setLevel(logging.WARNING)

        # Prepare the model parameters and forecast details
        model = "GR4JCN"
        params = (0.529, -3.396, 407.29, 1.072, 16.9, 0.947)

        forecast_duration = 3
        ts = get_local_testdata(
            "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
        )
        rvc = get_local_testdata("gr4j_cemaneige/solution.rvc")

        # Make the hindcasts for each initialization date. Here we will extract
        # ESP forecasts for a given calendar date for the years in "included_years"
        # as hindcast dates. Each ESP hindcast uses all available data in the ts dataset,
        # so in this case we will have 56/57 members for each hindcast initialization
        # depending on the date that we start on. The "hindcasts" dataset contains
        # all of the flow data from the ESP hindcasts for the initialization dates.
        # The "qobs" dataset contains all qobs in the timeseries: Climpred will
        # sort it all out during its processing. Note that the format of these datasets
        # is tailor-made to be used in climpred, and thus has specific dimension names.

        hindcasts, qobs = make_ESP_hindcast_dataset(
            model_name=model,
            forecast_date=dt.datetime(1955, 6, 30),
            included_years=list(range(1957, 1959)),
            forecast_duration=forecast_duration,
            ts=ts,
            hrus=(
                GR4JCN.LandHRU(
                    area=4250.6, elevation=843.0, latitude=54.4848, longitude=-123.3659
                ),
            ),
            params=params,
            rvc=str(rvc),
        )

        # Once we have the correctly formatted datasets, Make the hindcast object for climpred
        hindcast_object = make_climpred_hindcast_object(hindcasts, qobs)

        # This function is used to convert to binary to see if yes/no forecast is larger than obs
        def pos(x):
            return x > 0  # Check for binary outcome

        # These three functions respectively compute the rank histogram,
        # the crps and the reliability for the set of initialized dates
        # (i.e. forecast issue dates, here 1 day per year at the same calendar day).
        rank_histo_verif = hindcast_object.verify(
            metric="rank_histogram",
            comparison="m2o",
            dim=["member", "init"],
            alignment="same_inits",
        )
        crps_verif = hindcast_object.verify(
            metric="crps",
            comparison="m2o",
            dim=["member", "init"],
            alignment="same_inits",
        )
        reliability_verif = hindcast_object.verify(
            metric="reliability",
            comparison="m2o",
            dim=["member", "init"],
            alignment="same_inits",
            logical=pos,
        )

        assert "flow" in rank_histo_verif
        assert "flow" in crps_verif
        assert "flow" in reliability_verif
        assert rank_histo_verif.flow.shape[0] == forecast_duration
        assert reliability_verif.flow.shape[0] == forecast_duration
        assert crps_verif.flow.shape[0] == forecast_duration
