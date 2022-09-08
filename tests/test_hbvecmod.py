from ravenpy.extractors.routing_product import RoutingProductShapefileExtractor
from ravenpy.models import get_average_annual_runoff
from ravenpy.models.emulators.hbvecmod import HBVECMOD, HBVECMOD_OST
from ravenpy.utilities.testdata import get_local_testdata

TS = get_local_testdata("famine/famine_input.nc")
area = 100

default = [
    1.0,
    1.0,
    0.21941,
    0.15725,
    2.65,
    0.0,
    1.0,
    4.0,
    0.0464,
    1.0,
    1.0,
    1.0,
    1.0,
    0.01,
    0.01,
    1.0,
    0.03,
    0.03,
    1.1,
    0.02,
    100.0,
    5,
    1,
    0.1,
    1.0,
    0.1,
    0.01,
]

extractor = RoutingProductShapefileExtractor(
    get_local_testdata("famine/hru_Famine_final.zip"),
    routing_product_version="2.1",
)

bnds = [
    [1, 4],  # X1
    [1, 2],
    [0.05, 0.25],
    [0.05, 0.75],
    [0.05, 100],  # X5
    [-1, 1],
    [0, 4],
    [0, 7],
    [0.04, 0.07],
    [0, 5],  # X10
    [0.5, 2],
    [1, 7],
    [0.5, 2],
    [0.01, 1],
    [0.01, 1],  # X15
    [0.01, 10],
    [0.01, 1],
    [0.05, 0.1],
    [0.5, 2],
    [0.02, 0.2],  # X20
    [0.01, 100],
    [0, 5],
    [0.01, 1],
    [0.001, 1],
    [1, 3],  # X25
    [0.001, 1],
    [0.001, 5],
]


def test_simple():
    model = HBVECMOD()
    model.config.rvp.params = HBVECMOD.Params(*default)
    model.config.rvp.avg_annual_runoff = get_average_annual_runoff(
        TS, area * 1e6, obs_var="qobs"
    )
    rv_objs = extractor.extract()
    model.config.rvp.channel_profiles = rv_objs.pop("channel_profiles")

    for k, v in rv_objs.items():
        model.config.rvh.update(k, v)

    for sb in model.config.rvh.subbasins:
        sb.gauged = sb.subbasin_id == 160

    model.config.rvt.nc_index = [0, 1, 2]
    model.config.rvt.hydro_idx = (1,)
    model(ts=[TS], overwrite=True)


def test_calib_simple():
    model = HBVECMOD_OST()
    low, high = zip(*bnds)
    model.config.ost.lowerBounds = HBVECMOD.Params(*low)
    model.config.ost.upperBounds = HBVECMOD.Params(*high)
    model.config.rvp.avg_annual_runoff = get_average_annual_runoff(
        TS, area * 1e6, obs_var="qobs"
    )
    rv_objs = extractor.extract()
    model.config.rvp.channel_profiles = rv_objs.pop("channel_profiles")

    for k, v in rv_objs.items():
        model.config.rvh.update(k, v)

    for sb in model.config.rvh.subbasins:
        sb.gauged = sb.subbasin_id == 160

    model.config.rvt.nc_index = [0, 1, 2]
    model.config.rvt.hydro_idx = (1,)

    model.config.ost.max_iterations = 10
    model(ts=[TS], overwrite=True)
    assert model.calibrated_params
