from pathlib import Path

import pytest
import xarray

from ravenpy.utilities.testdata import (
    _default_cache_dir,
    get_file,
    open_dataset,
    query_folder,
)


class TestRemoteFileAccess:
    dap_url = "http://test.opendap.org:80/opendap/data/nc/"
    git_url = "https://github.com/Ouranosinc/raven-testdata"
    branch = "master"

    def test_get_file(self):
        file = get_file(name="ostrich-hbv-ec/raven-hbv-salmon.rvi", branch=self.branch)

        assert Path(_default_cache_dir).exists()
        assert file.is_file()
        with file.open() as f:
            header = f.read()
            assert ":FileType          rvi ASCII Raven 2.8.2" in header

    def test_open_dataset(self):
        ds = open_dataset(
            name="raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc",
            branch=self.branch,
        )

        assert (
            Path(_default_cache_dir)
            .joinpath(
                self.branch,
                "raven-gr4j-cemaneige",
                "Salmon-River-Near-Prince-George_meteo_daily.nc",
            )
            .exists()
        )
        assert isinstance(ds, xarray.Dataset)

    def test_open_dataset_no_cache(self):
        ds = open_dataset(
            name="raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
            branch=self.branch,
            cache=False,
        )

        assert (
            not Path(_default_cache_dir)
            .joinpath(
                "raven-gr4j-cemaneige",
                "Salmon-River-Near-Prince-George_meteo_daily_3d.nc",
            )
            .exists()
        )
        assert isinstance(ds, xarray.Dataset)

    @pytest.mark.xfail(reason="test.opendap.org is offline")
    def test_dap_access(self):
        ds = open_dataset(
            name="20070917-MODIS_A-JPL-L2P-A2007260000000.L2_LAC_GHRSST-v01.nc",
            dap_url=self.dap_url,
        )

        assert isinstance(ds, xarray.Dataset)


class TestQueryFolder:
    git_url = "https://github.com/Ouranosinc/raven-testdata"
    branch = "master"

    # def test_get_files_with_testdata_folder(self):
    #     rvs = get_file(TESTDATA["raven-gr4j-cemaneige-nc-rv"], branch="master")
    #     assert len(rvs) == 5

    # def test_query_folder(self):
    #     all_files = query_folder()
    #     assert len(all_files) == 129

    def test_query_specific_folder(self):
        folder = query_folder(folder="raven-gr4j-cemaneige", branch=self.branch)
        assert len(folder) == 8

    def test_query_folder_patterns(self):
        mohyse = query_folder(
            folder="/regionalisation_data/tests/", pattern="MOHYSE", branch=self.branch
        )
        assert len(mohyse) == 1
        assert mohyse[0] == str(
            Path("regionalisation_data", "tests", "MOHYSE_parameters.csv")
        )

    def test_query_folder_patterns_excessive_slashes(self):
        mohyse = query_folder(
            folder="///regionalisation_data/////tests///",
            pattern="MOHYSE",
            branch=self.branch,
        )
        assert len(mohyse) == 1
        assert mohyse[0] == str(
            Path("regionalisation_data", "tests", "MOHYSE_parameters.csv")
        )
