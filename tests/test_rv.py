import datetime as dt
import re
from collections import namedtuple
from io import StringIO
from pathlib import Path

import geopandas
import pytest

import ravenpy
from ravenpy.models.commands import GriddedForcingCommand
from ravenpy.models.importers import (
    RoutingProductGridWeightImporter,
    RoutingProductShapefileImporter,
)
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
from ravenpy.utilities.testdata import get_local_testdata


class TestRVFile:
    def test_simple_rv(self):
        fn = get_local_testdata("raven-hmets/*.rvp")
        rvf = RVFile(fn)

        assert rvf.ext == "rvp"
        assert rvf.stem == "raven-hmets-salmon"
        assert not rvf.is_tpl

    def test_simple_tpl(self):
        fn = get_local_testdata("ostrich-gr4j-cemaneige/*.rvp.tpl")
        rvf = RVFile(fn)

        assert rvf.ext == "rvp"
        assert rvf.stem == "raven-gr4j-salmon"
        assert rvf.is_tpl

    def test_ostIn(self):
        fn = get_local_testdata("ostrich-gr4j-cemaneige/ostIn.txt")
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
            unit="degC",
            dimensions=["time"],
        )
        tmp = str(v)

        assert compare(
            tmp,
            """:Data TEMP_MIN degC
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
            unit="degC",
            dimensions=["time"],
            linear_transform=(24000.0, 0.0),
        )

        assert ":LinearTransform 24000.000000000000000 0.000000000000000" in str(v)

    def test_deaccumulate(self):
        v = RavenNcData(
            var="tasmin",
            path="/path/tasmin.nc",
            var_name="tn",
            unit="degC",
            dimensions=["time"],
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
        rvc = open(get_local_testdata("gr4j_cemaneige/solution.rvc")).read()
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
        shp = get_local_testdata("raven-routing-sample/finalcat_hru_info.zip")
        importer = RoutingProductShapefileImporter(shp)
        sbs, land_group, lake_group, reservoirs, _, hrus = importer.extract()
        self.rvh = RVH(sbs, land_group, lake_group, reservoirs, hrus)

    def test_import_process(self):
        assert len(self.rvh.subbasins) == 46
        assert len(self.rvh.land_subbasin_group) == 41
        assert len(self.rvh.lake_subbasin_group) == 5
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
        _, _, _, _, cps, _ = importer.extract()
        self.rvp = RVP(cps)

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
        gfc = GriddedForcingCommand(grid_weights=gws)
        self.rvt = RVT([gfc])

    def test_import_process(self):
        res = self.rvt.to_rv()

        assert ":NumberHRUs 51" in res
        assert ":NumberGridCells 100" in res
        assert len(res.split("\n")) == 225


def test_isinstance_namedtuple():
    X = namedtuple("params", "x1, x2, x3")
    x = X(1, 2, 3)
    assert isinstance_namedtuple(x)
    assert not isinstance_namedtuple([1, 2, 3])
