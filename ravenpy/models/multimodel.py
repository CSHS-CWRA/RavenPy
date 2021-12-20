from pathlib import Path
from typing import List

from .base import Raven
from .emulators import get_model


class RavenMultiModel(Raven):
    def __init__(self, models, workdir=None):
        """Create multi-model raven instance.

        Parameters
        ----------
        models : sequence
          Model identifiers ('gr4jcn', 'hmets', 'mohyse', 'hbvec').
        """
        import tempfile

        self._names = models
        self._models = []

        workdir = workdir or tempfile.mkdtemp()
        Raven.__init__(self, workdir)

        for name in models:
            m = get_model(name)(workdir)
            m.model_dir = name
            self._models.append(m)

    def _rename_run_name(self, run_name=None):
        rns = {m.config.rvi.run_name for m in self._models}
        if (run_name is not None) or (len(rns) < len(self._models)):
            for m in self._models:
                rn = run_name or m.config.rvi.run_name
                m.config.rvi.run_name = rn + "-" + m.identifier

    def resume(self, solution=None):
        # TODO: Add support for model dependent solutions.
        for m in self._models:
            m.resume(solution)

    def run(self, ts, overwrite=False, **kwds):
        """Run model.

        Parameters
        ----------
        kwds : dict
          model_name : array
            Parameter array.
        """
        if overwrite:
            self.setup(overwrite)

        self._rename_run_name(kwds.pop("run_name", None))

        p = {}
        for m in self._models:
            p[m.identifier] = kwds.pop(m.identifier, None)

        procs = []
        for m in self._models:
            # Add params to kwds if passed in run.
            kw = kwds.copy()
            if p[m.identifier]:
                kw["params"] = p[m.identifier]

            procs.extend(m.run(ts, **kw))

        return procs

    def parse_results(self):
        # The Raven parent class uses `run_name` as the glob prefix, but here
        # since we have multiple models (with each its own config.rvi.run_name)
        # we simply use no prefix, which has the effect of getting all the files,
        # for all models
        patterns = {
            "hydrograph": "*Hydrographs.nc",
            "storage": "*WatershedStorage.nc",
            "solution": "*solution.rvc",
            "diagnostics": "*Diagnostics.csv",
        }

        for key, pattern in patterns.items():
            # There are no diagnostics if a streamflow time series is not provided.
            try:
                fns = self._get_output(pattern, path=self.exec_path)
            except UserWarning as exc:
                if key != "diagnostics":
                    raise exc
                else:
                    continue

            fns.sort()
            self.ind_outputs[key] = fns
            self.outputs[key] = self._merge_output(fns, pattern[1:])

        rv_paths = []
        for m in self._models:
            rv_paths += m._rv_paths

        self.outputs["rv_config"] = self._merge_output(rv_paths, "rv.zip")
