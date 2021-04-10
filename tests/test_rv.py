import datetime as dt
import re
from collections import namedtuple
from pathlib import Path

import pytest

import ravenpy
from ravenpy.config.commands import (
    BaseValueCommand,
    GriddedForcingCommand,
    MonthlyAverageCommand,
    RainCorrectionCommand,
)
from ravenpy.config.importers import (
    RoutingProductGridWeightImporter,
    RoutingProductShapefileImporter,
)
from ravenpy.config.rvs import RVC, RVH, RVI, RVP, RVT, Ost
from ravenpy.utilities.testdata import get_local_testdata


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


class TestOst:
    def test_random(self):
        o = Ost(None)
        assert o.random_seed == ""

        o.random_seed = 0
        assert o.random_seed == "RandomSeed 0"


class TestRVI:
    def test_supress_output(self):
        rvi = RVI(None)
        rvi.suppress_output = True
        assert rvi.suppress_output == ":SuppressOutput\n:DontWriteWatershedStorage"

        rvi = RVI(None)
        rvi.suppress_output = False
        assert rvi.suppress_output == ""


class TestRVC:
    @classmethod
    def setup_class(self):
        sol = open(get_local_testdata("gr4j_cemaneige/solution.rvc")).read()
        self.rvc = RVC(None)
        self.rvc.set_from_solution(sol)

    def test_parse(self):
        assert self.rvc.hru_state.atmosphere == 821.98274
        assert self.rvc.basin_state.qout == (1, 13.2166, 13.29232)

    def test_write(self):
        res = self.rvc.to_rv()
        assert "1,0.0,821.98274" in res
        assert ":BasinIndex 1 watershed" in res


class TestRVH:
    @classmethod
    def setup_class(self):
        shp = get_local_testdata("raven-routing-sample/finalcat_hru_info.zip")
        importer = RoutingProductShapefileImporter(shp)
        config = importer.extract()
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
    @classmethod
    def setup_class(self):
        shp = get_local_testdata("raven-routing-sample/finalcat_hru_info.zip")
        importer = RoutingProductShapefileImporter(shp)
        config = importer.extract()
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
    @classmethod
    def setup_class(self):
        input_file = get_local_testdata("raven-routing-sample/VIC_streaminputs.nc")
        routing_file = get_local_testdata("raven-routing-sample/finalcat_hru_info.zip")
        importer = RoutingProductGridWeightImporter(input_file, routing_file)
        gws = importer.extract()
        self.gfc = GriddedForcingCommand(grid_weights=gws)

    def test_import_process(self):
        res = self.gfc.to_rv()

        assert ":NumberHRUs 51" in res
        assert ":NumberGridCells 100" in res
        # FIXME: This test is not superb.
        assert len(res.split("\n")) == 226


class TestBaseValueCommand:
    def test_raincorrection(self):
        rc = RainCorrectionCommand(3)
        assert f"{rc}" == ":RainCorrection 3"
