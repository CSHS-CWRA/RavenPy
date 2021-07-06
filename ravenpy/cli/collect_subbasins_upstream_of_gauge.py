from collections import defaultdict
from pathlib import Path

import click
import geopandas as gpd


@click.command()
@click.argument("input-file", type=click.Path(exists=True))
@click.argument("gauge-id", type=str)
@click.option("-o", "--output", type=click.Path(), help="Name of the output shapefile.")
def collect_subbasins_upstream_of_gauge(
    input_file,
    gauge_id,
    output,
):
    """
    Find the subbasins upstream of a gauge from a Routing Product shapefile, and
    save them in a new shapefile.

    INPUT_FILE: Routing Product shapefile (e.g. "drainage_region_0003_v2-1/finalcat_info_v2-1.shp").

    GAUGE_ID: ID of the target gauge, to be found in the "Obs_NM" column (e.g. "02LE024").
    """
    if Path(input_file).suffix == ".zip":
        input_file = f"zip://{input_file}"
    gdf = gpd.read_file(input_file)

    if gauge_id not in gdf.Obs_NM.values:
        raise click.ClickException(f"Cannot find gauge {gauge_id} in `Obs_NM` column")

    if not output:
        p = Path(input_file)
        output_file = p.parent / f"{p.stem}_upstream_of_gauge_{gauge_id}.shp"
    else:
        output_file = output

    downsubid_to_subids = defaultdict(set)

    for _, r in gdf.iterrows():
        downsubid_to_subids[r.DowSubId].add(r.SubId)

    # Starting from the SubId of the gauge, iteratively expand the set of upstream subbasin SubIds
    upstream_subids = {gdf[gdf.Obs_NM == gauge_id].iloc[0].SubId}
    prev = upstream_subids.copy()
    while True:
        curr = set()
        for did in prev:
            curr |= downsubid_to_subids[did]
        if curr:
            upstream_subids |= curr
            prev = curr
        else:
            break

    gdf_upstream = gdf[gdf.SubId.isin(upstream_subids)]

    gpd.GeoDataFrame(gdf_upstream).to_file(output_file)

    click.echo(
        f"Found {len(gdf_upstream)} upstream subbasins, saved them in {output_file}"
    )
