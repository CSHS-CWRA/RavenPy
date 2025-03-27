import tempfile
import warnings
from collections.abc import Sequence
from dataclasses import astuple
from pathlib import Path
from typing import Optional, Union

from spotpy.parameter import Uniform, generate

from ravenpy import Emulator
from ravenpy.config.base import Params
from ravenpy.config.options import evaluation_metrics_multiplier
from ravenpy.config.rvs import Config


class SpotSetup:
    def __init__(
        self,
        config: Config,
        low: Union[Params, Sequence],
        high: Union[Params, Sequence],
        workdir: Optional[Union[Path, str]] = None,
    ):
        """Class to configure spotpy with Raven emulators.

        Parameters
        ----------
        config : Config
            Emulator Config instance with symbolic expressions.
        low : Union[Params, Sequence]
            Lower boundary for parameters.
        high : Union[Params, Sequence]
            Upper boundary for parameters.
        workdir : Union[str, Path], optional
            Work directory. If None, a temporary directory will be created.
        """
        if not config.is_symbolic:
            raise ValueError(
                "config should be a symbolic configuration, where params are not set to their numerical "
                "values."
            )

        self.config = config
        self.path = Path(workdir or tempfile.mkdtemp())
        self.diagnostics = None

        if config.suppress_output is not True:
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
            pdist.append(  # noqa: PERF401
                Uniform(self.pnames[i], low=low[i], high=high[i])
            )
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

    @staticmethod
    def evaluation():
        """Return the observation.

        Since Raven computes the objective function itself, we simply return a placeholder.
        """
        return 1

    def simulation(self, x):
        """Run the model, but return a placeholder value instead of the model output."""
        self._iteration += 1

        # Update parameters
        c = self.config.set_params(list(x))

        # Create emulator instance
        emulator = Emulator(config=c, workdir=self.path / f"c{self._iteration:03}")

        # Run the model
        output = emulator.run()

        self.diagnostics = output.diagnostics

        self._iteration += 1

        return 1

    # FIXME: Verify that this function is correct
    def objectivefunction(self, evaluation, simulation):  # noqa: F841
        """Return the objective function.

        Note that we short-circuit the evaluation and simulation entries, since the objective function has already
        been computed by Raven.
        """
        out = [
            self.diagnostics[f"DIAG_{m}"][0] * self._multipliers[m]
            for m in self.metrics
        ]

        return out
