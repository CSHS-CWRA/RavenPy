import datetime as dt
import re

import pytest

from ravenpy.config.commands import EvaluationPeriod, GriddedForcingCommand
from ravenpy.config.rvs import OST, RVC, RVH, RVI, RVP, RVT, Config
from ravenpy.extractors import (
    RoutingProductGridWeightExtractor,
    RoutingProductShapefileExtractor,
)


class TestRV:
    def test_end_date(self):
        rvi = RVI(None)
        rvi.run_name = "test"
        rvi.start_date = dt.datetime(2000, 1, 1)
        rvi.end_date = dt.datetime(2000, 1, 11)

        assert 10 == rvi.duration

        rvi.duration = 11
        assert dt.datetime(2000, 1, 12) == rvi.end_date

    def test_evaluation_metrics(self):
        rvi = RVI(None)
        rvi.evaluation_metrics = "LOG_NASH"

        with pytest.raises(ValueError):
            rvi.evaluation_metrics = "JIM"

    def test_evaluation_periods(self):
        rvi = RVI(None)
        assert rvi.evaluation_periods == ""

        rvi.evaluation_periods = [
            EvaluationPeriod("dry", "1980-01-01", "1989-12-31"),
            EvaluationPeriod("wet", "1990-01-01", "2000-12-31"),
        ]
        out = rvi.evaluation_periods
        assert len(out.split("\n")) == 2
        assert out.startswith(":EvaluationPeriod")

        # Check date input
        d = EvaluationPeriod("dry", dt.date(1980, 1, 1), dt.date(1989, 12, 31))
        assert str(d) == str(rvi.evaluation_periods.splitlines()[0])


class TestRVI:
    def test_supress_output(self):
        rvi = RVI(None)
        rvi.suppress_output = True
        assert rvi.suppress_output == ":SuppressOutput\n:DontWriteWatershedStorage\n"

        rvi = RVI(None)
        rvi.suppress_output = False
        assert rvi.suppress_output == ""


class TestRVC:
    @pytest.fixture(autouse=True)
    def setup(self, get_file):
        sol = open(get_file("gr4j_cemaneige/solution.rvc")).read()
        self.rvc = RVC.create_solution(sol)

    def test_parse(self):
        assert len(self.rvc.hru_states) == 1
        assert self.rvc.hru_states[1].data["ATMOSPHERE"] == 821.98274
        assert self.rvc.hru_states[1].data["ATMOS_PRECIP"] == -1233.16

        assert len(self.rvc.basin_states) == 1
        assert self.rvc.basin_states[1].channel_storage == 0
        assert self.rvc.basin_states[1].qout == (1, 13.21660, 13.29232)

    def test_format(self):
        rv = self.rvc.to_rv()
        assert ":BasinIndex 1 watershed" in rv


class TestRVH:
    @pytest.fixture(autouse=True)
    def setup(self, get_file):
        shp = get_file("raven-routing-sample/finalcat_hru_info.zip")
        extractor = RoutingProductShapefileExtractor(shp)
        config = extractor.extract()
        self.rvh = RVH(None)
        for k, v in config.items():
            if k != "channel_profiles":
                self.rvh.update(k, v)

    def test_import_process(self):
        assert len(self.rvh.subbasins) == 46
        assert len(self.rvh.land_subbasin_ids) == 41
        assert len(self.rvh.lake_subbasin_ids) == 5
        assert len(self.rvh.reservoirs) == 5
        assert len(self.rvh.hrus) == 51

    def test_format(self):
        res = self.rvh.to_rv()

        sbs = (
            re.search(":SubBasins(.+):EndSubBasins", res, re.MULTILINE | re.DOTALL)
            .group(1)
            .split("\n")
        )
        sbs = list(filter(None, sbs))  # remove whitespaces
        assert len(sbs) == len(self.rvh.subbasins) + 2

        assert res.count("ZERO-") == len(self.rvh.reservoirs)

        hrus = (
            re.search(":HRUs(.+):EndHRUs", res, re.MULTILINE | re.DOTALL)
            .group(1)
            .split("\n")
        )
        hrus = list(filter(None, hrus))  # remove whitespaces
        assert len(hrus) == len(self.rvh.hrus) + 2

        assert res.count(":Reservoir") == len(self.rvh.reservoirs)


class TestRVP:
    @pytest.fixture(autouse=True)
    def setup(self, get_file):
        shp = get_file("raven-routing-sample/finalcat_hru_info.zip")
        extractor = RoutingProductShapefileExtractor(shp)
        config = extractor.extract()
        self.rvp = RVP(None)
        self.rvp.tmpl = "{channel_profiles}"
        self.rvp.channel_profiles = config["channel_profiles"]

    def test_import_process(self):
        assert len(self.rvp.channel_profiles) == 46

    def test_format(self):
        res = self.rvp.to_rv()

        assert res.count(":ChannelProfile") == 46
        assert res.count(":EndChannelProfile") == 46


class TestRVT:
    @pytest.fixture(autouse=True)
    def setup(self, get_file):
        input_file = get_file("raven-routing-sample/VIC_streaminputs.nc")
        routing_file = get_file("raven-routing-sample/finalcat_hru_info.zip")
        extractor = RoutingProductGridWeightExtractor(input_file, routing_file)
        gws = extractor.extract()
        self.gfc = GriddedForcingCommand(grid_weights=gws)

    def test_import_process(self):
        res = self.gfc.to_rv()

        assert ":NumberHRUs 51" in res
        assert ":NumberGridCells 100" in res
        # FIXME: This test is not superb.
        assert len(res.split("\n")) == 224

    def test_gauges(self, get_file):
        rvt = RVT(config=None)
        rvt.meteo_idx = [0, 1, 2]
        rvt.configure_from_nc_data([get_file("famine/famine_input.nc")])
        out = rvt.to_rv()
        assert len(re.findall(":Gauge", out)) == 3
        assert set(re.findall(r":StationIdx\s*(\d)", out)) == {"1", "2", "3"}