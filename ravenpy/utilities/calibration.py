#!/usr/bin/env python3
"""
Created on Sat Apr 30 16:18:54 2022

@author: ets
"""

import warnings
from dataclasses import astuple
from pathlib import Path
from typing import Sequence, Union

import numpy as np
from spotpy.parameter import Uniform, generate

from ravenpy import Emulator
from ravenpy.new_config.base import Params
from ravenpy.new_config.options import evaluation_metrics_multiplier
from ravenpy.new_config.rvs import Config
from ravenpy.ravenpy import parse_diagnostics


class SpotSetup:
    def __init__(
        self,
        config: Config,
        low: Union[Params, Sequence],
        high: [Params, Sequence],
        path: Union[Path, str],
    ):
        """
        Class to configure spotpy with Raven emulators.

        Parameters
        ----------
        config: Config
          A Config instance with symbolic expressions.
        """
        self.config = config
        self.path = path
        self.diagnostics = None

        if config.suppress_output is not False:
            warnings.warn(
                "Add the `SuppressOutput` command to the configuration to reduce IO."
            )

        if config.evaluation_metrics is None:
            raise AttributeError(":EvaluationMetrics is undefined.")

        # Get evaluation metrics and their multiplier (spotpy maximizes the obj function)
        self.metrics = [m.value for m in config.evaluation_metrics]
        self._multipliers = {m: evaluation_metrics_multiplier[m] for m in self.metrics}

        p = config.params

        self.pnames = list(p.__dataclass_fields__.keys())
        self.pdist = self.init_params(low, high)
        self._iteration = 0

    def init_params(self, low: Union[Params, Sequence], high: [Params, Sequence]):
        # Validate parameters
        low = astuple(self._to_dataclass(low))
        high = astuple(self._to_dataclass(high))

        pdist = []
        for i in range(len(low)):
            pdist.append(Uniform(self.pnames[i], low=low[i], high=high[i]))
        return pdist

    def _to_dataclass(self, p):
        """Validate parameters against Params dataclass from model configuration."""
        kls = self.config.params.__class__
        if isinstance(p, kls):
            return p
        return kls(*p)

    def parameters(self):
        """Return a random parameter combination."""
        return generate(self.pdist)

    def evaluation(self):
        """Return the observation. Since Raven computes the objective function itself,
        we simply return a placeholder."""
        return 1

    def simulation(self, x):
        """Run the model, but return a placeholder value instead of the model output."""
        self._iteration += 1

        # Update parameters
        c = self.config.set_params(list(x))

        # Create emulator instance
        emulator = Emulator(config=c, path=self.path / f"c{self._iteration:03}")

        # Run the model
        emulator.build()
        emulator.run()

        self.diagnostics = emulator.parse_diagnostics()

        self._iteration += 1

        return 1

    def objectivefunction(self, evaluation, simulation):
        """Return the objective function.

        Note that we short-circuit the evaluation and simulation entries, since the objective function has already
        been computed by Raven.
        """

        return [
            self.diagnostics[f"DIAG_{m}"][0] * self._multipliers[m]
            for m in self.metrics
        ]


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
