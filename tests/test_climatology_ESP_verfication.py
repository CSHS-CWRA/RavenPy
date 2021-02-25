#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 23 22:25:32 2021

@author: ets
"""

import datetime as dt

import matplotlib.pyplot as plt


from ravenpy.utilities.testdata import get_local_testdata
from ravenpy.utilities.forecasting import make_ESP_hindcast_dataset


class TestClimatologyESP:
    def test_simple(self):
        
        model = "GR4JCN"
        params = (0.529, -3.396, 407.29, 1.072, 16.9, 0.947)

        ts=get_local_testdata("raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc")
        
        hindcasts,qobs=make_ESP_hindcast_dataset(model_name=model,
                                                              forecast_date=dt.datetime(1954, 12, 30),                                                              
                                                              included_years=list(range(1954,1956)),
                                                              forecast_duration=90,
                                                              ts=ts,  
                                                              area="4250.6",
                                                              elevation="843.0",
                                                              latitude=54.4848,
                                                              longitude=-123.3659,
                                                              params=params,
                                                              )
        

       
