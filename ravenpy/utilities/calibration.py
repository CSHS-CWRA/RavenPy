#!/usr/bin/env python3
"""
Created on Sat Apr 30 16:18:54 2022

@author: ets
"""

import numpy as np
from spotpy.parameter import Uniform, generate


class spotpy_setup:
    def __init__(self, model, ts, obj_func=None):
        # Just a way to keep this example flexible and applicable to various examples
        self.obj_func = obj_func

        # self.params=[Uniform('X1',low=1, high=10),
        #             Uniform('X2',low=1, high=10),
        #             Uniform('X3',low=1, high=10),
        #             Uniform('X4',low=1, high=10),
        #             Uniform('X5',low=1, high=10),
        #             Uniform('X6',low=1, high=10)
        #             ]
        self.params = []
        for i in range(0, len(model.low)):
            self.params.append(Uniform(str(i), low=model.low[i], high=model.high[i]))
        self.model = model
        self.ts = ts
        return

    def evaluation(self):
        return 1

    def parameters(self):
        return generate(self.params)

    def simulation(self, vector):
        x = np.array(vector)
        self.x = x
        return 1

    def model(self):
        return self

    def objectivefunction(self, evaluation, simulation):
        model = self.model

        # TODO: Add warm-up period!
        model.config.update("params", self.x)
        model(self.ts)
        d = model.diagnostics
        objfun = d["DIAG_NASH_SUTCLIFFE"][0]

        return objfun
