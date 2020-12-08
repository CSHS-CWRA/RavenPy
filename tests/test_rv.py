import datetime as dt
import re
from collections import namedtuple
from io import StringIO
from pathlib import Path

import geopandas
import pytest

import ravenpy
from ravenpy.models.importers import RoutingProductShapefileImporter
from ravenpy.models.rv import (
    RV,
    RVC,
    RVH,
    RVI,
    RVP,
    RVT,
    MonthlyAverage,
    Ost,
    RavenNcData,
    RVFile,
    isinstance_namedtuple,
)

from .common import TESTDATA


class TestRVFile:
    def test_simple_rv(self):
        fn = list(TESTDATA["raven-hmets"].glob("*.rvp"))[0]
        rvf = RVFile(fn)

        assert rvf.ext == "rvp"
        assert rvf.stem == "raven-hmets-salmon"
        assert not rvf.is_tpl

    def test_simple_tpl(self):
        fn = list(TESTDATA["ostrich-gr4j-cemaneige"].glob("*.rvp.tpl"))[0]
        rvf = RVFile(fn)

        assert rvf.ext == "rvp"
        assert rvf.stem == "raven-gr4j-salmon"
        assert rvf.is_tpl

    def test_ostIn(self):
        fn = list(TESTDATA["ostrich-gr4j-cemaneige"].glob("ostIn.txt"))[0]
        rvf = RVFile(fn)

        assert rvf.ext == "txt"
        assert rvf.stem == "ostIn"
        assert rvf.is_tpl

    def test_tags(self):
        rvp = list(
            (Path(ravenpy.__file__).parent / "models" / "raven-gr4j-cemaneige").glob(
                "*.rvp"
            )
        )[0]
        rvf = RVFile(rvp)

        assert isinstance(rvf.tags, list)
        assert "params.GR4J_X3" in rvf.tags

    def test_fail(self):
        fn = Path(ravenpy.__file__).parent
        with pytest.raises(ValueError):
            RVFile(fn)


class TestRV:
    def test_end_date(self):
        rvi = RVI(
            run_name="test",
            start_date=dt.datetime(2000, 1, 1),
            end_date=dt.datetime(2000, 1, 11),
        )

        assert 10 == rvi.duration

        rvi.duration = 11
        assert dt.datetime(2000, 1, 12) == rvi.end_date

    def test_params(self):
        class RVP(RV):
            params = namedtuple("p", "x, y")

        rvp = RVP()
        rvp.params = RVP.params(1, 2)
        assert rvp.params.x == 1

    def test_dict_interface(self):
        rv = RV(run_name="test")

        assert rv["run_name"] == rv.run_name

        with pytest.raises(AttributeError):
            rv["r"] = 6

    def test_evaluation_metrics(self):
        rvi = RVI()
        rvi.evaluation_metrics = "LOG_NASH"

        with pytest.raises(ValueError):
            rvi.evaluation_metrics = "JIM"

    def test_update(self):
        rv = RV(a=None, b=None)
        rv.update({"a": 1, "b": 2})
        assert rv.a == 1

        rv.c = 1
        assert rv["c"] == 1

    def test_namedtuple(self):
        class Mod(RV):
            params = namedtuple("params", "x1, x2, x3")

        m = Mod(params=Mod.params(1, 2, 3))
        assert m.params.x1 == 1


def compare(a, b):
    """
    Compare two base strings, disregarding whitespace
    """
    import re

    return re.sub(r"\s*", "", a) == re.sub(r"\s*", "", b)


class TestRavenNcData:
    def test_simple(self):
        v = RavenNcData(
            var="tasmin",
            path="/path/tasmin.nc",
            var_name="tn",
            unit="deg_C",
            dimensions=[
                "time",
            ],
        )
        tmp = str(v)

        assert compare(
            tmp,
            """:Data TEMP_MIN deg_C
                                  :ReadFromNetCDF
                                     :FileNameNC      /path/tasmin.nc
                                     :VarNameNC       tn
                                     :DimNamesNC      time
                                     :StationIdx      1
                                  :EndReadFromNetCDF
                               :EndData""",
        )

    def test_linear_transform(self):
        v = RavenNcData(
            var="tasmin",
            path="/path/tasmin.nc",
            var_name="tn",
            unit="deg_C",
            dimensions=[
                "time",
            ],
            linear_transform=(24000.0, 0.0),
        )

        assert ":LinearTransform 24000.000000000000000 0.000000000000000" in str(v)

    def test_deaccumulate(self):
        v = RavenNcData(
            var="tasmin",
            path="/path/tasmin.nc",
            var_name="tn",
            unit="deg_C",
            dimensions=[
                "time",
            ],
            deaccumulate=True,
        )

        assert ":Deaccumulate" in str(v)


class TestMonthlyAve:
    def test_simple(self):
        ave = str(MonthlyAverage("Evaporation", range(12)))
        assert ave.startswith(":MonthlyAveEvaporation, 0, 1, 2")


class TestOst:
    def test_random(self):
        o = Ost()
        assert o.random_seed == ""

        o.random_seed = 0
        assert o.random_seed == "RandomSeed 0"


class TestRVI:
    def test_supress_output(self):
        rvi = RVI(suppress_output=True)
        assert rvi.suppress_output == ":SuppressOutput\n:DontWriteWatershedStorage"

        rvi = RVI(suppress_output=False)
        assert rvi.suppress_output == ""


class TestRVC:
    @classmethod
    def setup_class(self):
        rvc = TESTDATA["solution.rvc"].read_text()
        self.r = RVC()
        self.r.parse(rvc)

    def test_parse(self):
        assert self.r.hru_state.atmosphere == 821.98274
        assert self.r.basin_state.qout == [
            13.21660,
        ]
        assert self.r.basin_state.qoutlast == 13.29232

    def test_write(self):
        assert self.r.txt_hru_state.startswith("1,")
        assert self.r.txt_basin_state.strip().startswith(":BasinIndex 1,watershed")

    def test_format(self):
        rvc_template = Path(ravenpy.models.__file__).parent / "global" / "global.rvc"
        params = dict(self.r.items())
        rvc_template.read_text().format(**params)


class TestRVH:
    @classmethod
    def setup_class(self):
        importer = RoutingProductShapefileImporter(
            f"zip://{TESTDATA['routing-sample']}"
        )
        sbs, groups, lakes, _, hrus = importer.extract()
        self.rvh = RVH(sbs, groups, lakes, hrus)

    def test_import_process(self):
        assert len(self.rvh._subbasins) == 46
        assert len(self.rvh._land_subbasin_group_ids) == 41
        assert len(self.rvh._lake_subbasin_group_ids) == 5
        assert len(self.rvh._lakes) == 5
        assert len(self.rvh._hrus) == 51

    def test_format(self):
        rvh_template = Path(ravenpy.models.__file__).parent / "global" / "global.rvh"
        params = dict(self.rvh.items())
        res = rvh_template.read_text().format(**params)

        sbs = (
            re.search(":SubBasins(.+):EndSubBasins", res, re.MULTILINE | re.DOTALL)
            .group(1)
            .split("\n")
        )
        sbs = list(filter(None, sbs))  # remove whitespaces
        assert len(sbs) == len(self.rvh._subbasins) + 2

        assert res.count("ZERO-") == len(self.rvh._lakes)

        hrus = (
            re.search(":HRUs(.+):EndHRUs", res, re.MULTILINE | re.DOTALL)
            .group(1)
            .split("\n")
        )
        hrus = list(filter(None, hrus))  # remove whitespaces
        assert len(hrus) == len(self.rvh._hrus) + 2

        assert res.count(":Reservoir") == len(self.rvh._lakes)


class TestRVP:
    @classmethod
    def setup_class(self):
        importer = RoutingProductShapefileImporter(
            f"zip://{TESTDATA['routing-sample']}"
        )
        _, _, _, cps, _ = importer.extract()
        self.rvp = RVP(cps)

    def test_import_process(self):
        assert len(self.rvp._channel_profiles) == 46

    def test_format(self):
        res = self.rvp.render_to_rv()

        assert res.count(":ChannelProfile") == 46
        assert res.count(":EndChannelProfile") == 46


def test_isinstance_namedtuple():
    X = namedtuple("params", "x1, x2, x3")
    x = X(1, 2, 3)
    assert isinstance_namedtuple(x)
    assert not isinstance_namedtuple([1, 2, 3])
