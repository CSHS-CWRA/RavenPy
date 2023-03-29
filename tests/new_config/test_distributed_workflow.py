from ravenpy.extractors.new_config.routing_product import (
    BasinMakerExtractor,
    open_shapefile,
    upstream_from_coords,
)
from ravenpy.new_config import commands as rc

"""
Distributed model workflow

"""


def test_simple_workflow(get_local_testdata, minimal_emulator, tmp_path):
    shp_path = get_local_testdata(
        "basinmaker/drainage_region_0175_v2-1/finalcat_info_v2-1.zip"
    )

    # Note that for this to work, the coordinates must be in the small
    # BasinMaker example (drainage_region_0175)
    df = open_shapefile(shp_path)

    # Gauge station for observations at Matapedia
    # SubId: 175000128
    # -67.12542 48.10417
    sub = upstream_from_coords(-67.12542, 48.10417, df)

    # Confirm we got the right watershed
    assert 175000128 in sub["SubId"].to_list()

    bm = BasinMakerExtractor(
        df=sub,
        hru_aspect_convention="ArcGIS",
    )
    rvh = bm.extract()

    obs = rc.ObservationData.from_nc(
        get_local_testdata("matapedia/Qobs_Matapedia_01BD009.nc"),
        alt_names=("discharge",),
    )

    conf = minimal_emulator.copy(update={**rvh, "ObservationData": obs})
    conf.write_rv(workdir=tmp_path)
