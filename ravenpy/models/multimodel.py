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
                m.config.rvi.run_name = rn + "_" + m.config.identifier

    def assign(self, key, value):
        """Assign key to all models, unless it's model parameters."""
        # Model parameter case
        if key in self._names:
            m = self._models[self._names.index(key)]
            m.assign("params", value)
        else:
            for m in self._models:
                m.assign(key, value)

    def resume(self, solution=None):
        # TODO: Add support for model dependent solutions.
        for m in self._models:
            m.resume(solution)

    @property
    def _rv_paths(self):
        out = []
        for m in self._models:
            out.extend(m._rv_paths)
        return out

    @_rv_paths.setter
    def _rv_paths(self, value):
        pass

    # @property
    # def rvs(self):
    #     out = []
    #     for m in self._models:
    #         out.extend(m._rv_paths)
    #     return out

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
            p[m.config.identifier] = kwds.pop(m.config.identifier, None)

        procs = []
        for m in self._models:
            # Add params to kwds if passed in run.
            kw = kwds.copy()
            if p[m.config.identifier]:
                kw["params"] = p[m.config.identifier]

            procs.extend(m.run(ts, **kw))

        return procs
