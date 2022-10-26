#!/usr/bin/env python3
"""
Created on Sat Apr 30 16:18:54 2022

@author: ets
"""

import numpy as np
from spotpy.parameter import Uniform, generate


class SpotpySetup:
    def __init__(self, model, ts, obj_func=None):
        """

        Parameters
        ----------
        model: Raven emulator subclass instance
          Raven emulator.
        ts: list
          Forcing files.
        obj_func: func
          Objective function.
        """
        self.model = model

        # Make sure no output is written to disk
        self.model.config.rvi.suppress_output = True

        self.ts = ts

        # Just a way to keep this example flexible and applicable to various examples
        self.obj_func = obj_func

        # Initialize parameters
        self.params = []
        for i in range(0, len(model.low)):
            self.params.append(Uniform(str(i), low=model.low[i], high=model.high[i]))

    def evaluation(self):
        """In theory this method should return the true value. Since Raven computes the objective function,
        we simply return a placeholder."""
        return 1

    def parameters(self):
        """Return a random parameter combination."""
        return generate(self.params)

    def simulation(self, x):
        """Run the model, but return a placeholder value."""

        # Update parameters
        self.model.config.update("params", np.array(x))

        # Run the model
        self.model._execute(self.ts)
        return 1

    def objectivefunction(self, evaluation, simulation, params=None):
        """Return the objective function.

        Note that we short-circuit the evaluation and simulation entries, since the objective function has already
        been computed by Raven.
        """
        d = self.model.diagnostics
        return d["DIAG_NASH_SUTCLIFFE"][0]
