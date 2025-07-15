"""Distributed model workflow."""

import datetime as dt

from ravenpy import Emulator
from ravenpy.config import commands as rc
from ravenpy.config.emulators import GR4JCN
from ravenpy.extractors.routing_product import (
    BasinMakerExtractor,
    GridWeightExtractor,
    open_shapefile,
    upstream_from_coords,
)


def test_simple_workflow(tmp_path, yangtze):
    shp_path = yangtze.fetch(
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

    # Extract the subbasins and HRUs (one HRU per sub-basin)
    bm = BasinMakerExtractor(
        df=sub,
        hru_aspect_convention="ArcGIS",
    )
    rvh = bm.extract(hru_from_sb=True)

    # Streamflow obs
    qobs_fn = yangtze.fetch("matapedia/Qobs_Matapedia_01BD009.nc")

    qobs = rc.ObservationData.from_nc(
        qobs_fn,
        uid=175000128,
        alt_names=("discharge",),
    )

    # Meteo obs for GriddedForcing - does not work because subbasins do not overlap 100% with the ERA data
    meteo_grid_fn = yangtze.fetch("matapedia/Matapedia_meteo_data_2D.nc")

    # Dict of GW attributes
    gw = GridWeightExtractor(
        meteo_grid_fn,
        shp_path,
        dim_names=("longitude", "latitude"),
        var_names=("longitude", "latitude"),
        gauge_ids=[
            "01BD009",
        ],
    ).extract()

    assert gw["number_hrus"] == len(sub)

    # Write GW command to file
    gw_fn = tmp_path / "gw.txt"
    gw_fn.write_text(rc.GridWeights(**gw).to_rv())

    forcing = {"TEMP_MIN": "tmin", "TEMP_MAX": "tmax", "PRECIP": "pr"}

    [
        rc.GriddedForcing.from_nc(
            meteo_grid_fn, dtyp, alt_names=(alias,), grid_weights=gw_fn
        )
        for (dtyp, alias) in forcing.items()
    ]
    # Weights for some HRUs do not sum to one.

    # Meteo forcing per station (virtual stations, since this is ERA5 data)
    meteo_station = yangtze.fetch("matapedia/Matapedia_meteo_data_stations.nc")

    [
        rc.StationForcing.from_nc(meteo_station, dtyp, alt_names=(alias,))
        for (dtyp, alias) in forcing.items()
    ]
    # TODO: Complete with weights calculations

    # Virtual Gauges
    gauges = [
        rc.Gauge.from_nc(
            meteo_station,
            data_type=[s for s in forcing.keys()],
            station_idx=i + 1,
            alt_names=forcing,
        )
        for i in range(6)
    ]

    conf = GR4JCN(
        params=[0.529, -3.396, 407.29, 1.072, 16.9, 0.947],
        StartDate=dt.datetime(2000, 1, 1),
        Duration=15,
        GlobalParameter={"AVG_ANNUAL_RUNOFF": 208.480},
        ObservationData=[qobs],
        Gauge=gauges,
        **rvh,
    )

    out = Emulator(conf, workdir=tmp_path).run()

    # Number of gauged sub-basins
    ng = sum([sb.gauged for sb in conf.sub_basins])
    assert len(out.hydrograph.nbasins) == ng
