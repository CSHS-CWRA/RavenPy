import pytest
import xarray as xr

from ravenpy.config.commands import (
    LandUseClassesCommand,
    SBGroupPropertyMultiplierCommand,
    SoilClassesCommand,
    SoilProfilesCommand,
    VegetationClassesCommand,
)
from ravenpy.extractors.routing_product import (
    RoutingProductGridWeightExtractor,
    RoutingProductShapefileExtractor,
)
from ravenpy.models import Raven
from ravenpy.utilities.testdata import get_local_testdata


class TestRouting:
    def test_lievre_tutorial(self):
        """
        This test reproduces the Lievre tutorial setup:

        http://raven.uwaterloo.ca/files/RavenTutorial6.zip

        """

        ###############
        # Input files #
        ###############

        routing_product_shp_path = get_local_testdata(
            "raven-routing-sample/finalcat_hru_info.zip"
        )

        vic_streaminputs_nc_path = get_local_testdata(
            "raven-routing-sample/VIC_streaminputs.nc"
        )
        vic_temperatures_nc_path = get_local_testdata(
            "raven-routing-sample/VIC_temperatures.nc"
        )

        observation_data_nc_path = get_local_testdata(
            "raven-routing-sample/WSC02LE024.nc"
        )

        #########
        # Model #
        #########

        model = Raven()

        model.config.identifier = "raven-lievre-routing"

        #######
        # RVI #
        #######

        model.config.rvi.tmpl = """
        :CatchmentRoute        ROUTE_DUMP                    # Catchment routing method, used to convey water from the catchment tributaries and rivulets to the subbasin outlets. DEFAULT ROUTE_DUMP, which instantly ‘dumps’ all water in the subbasin stream reach.
        :Routing               ROUTE_DIFFUSIVE_WAVE          # Channel routing method which is used to transport water from upstream to downstream within the main subbasin channels. DEFAULT ROUTE_DIFFUSIVE_WAVE, which analytically solves the diffusive wave equation along the reach using a constant reference celerity.
        :PrecipIceptFract      PRECIP_ICEPT_NONE             # Estimation of the precipitation interception fraction. In this routing model, stream input(s) are "pretending" to be precipitation going into Raven, thus using DEFAULT PRECIP_ICEPT_NONE to indicate no interception processes are adopted.
        :PotentialMeltMethod   POTMELT_NONE                  # Estimation of the potential snow melt. In this routing model, snow melt processes are not relevant, thus using DEFAULT POTMELT_NONE method.
        :SoilModel             SOIL_ONE_LAYER                # In this routing model, use DEFAULT SOIL_ONE_LAYER to define single soil layer structure.

        :HydrologicProcesses
          :Precipitation     PRECIP_RAVEN             ATMOS_PRECIP     PONDED_WATER          # Moves stream input(s) from ATMOS_PRECIP to PONDED_WATER storage (waiting for runoff). Use DEFAULT PRECIP_RAVEN method.
          :Flush             RAVEN_DEFAULT            PONDED_WATER     SURFACE_WATER         # Moves water from PONDED_WATER to SURFACE_WATER (routed to outlet). Use DEFAULT RAVEN_DEFAULT method.
        :EndHydrologicProcesses
        """

        streaminputs = xr.open_dataset(vic_streaminputs_nc_path)

        start = streaminputs.indexes["time"][0]
        end = streaminputs.indexes["time"][-4]  # to match the tutorial end date

        model.config.rvi.start_date = start.to_pydatetime()
        model.config.rvi.end_date = end.to_pydatetime()

        # Raven will use 24h even though the NC inputs are 6h
        model.config.rvi.time_step = "24:00:00"

        model.config.rvi.evaluation_metrics = [
            "NASH_SUTCLIFFE",
            "PCT_BIAS",
            "KLING_GUPTA",
        ]

        #######
        # RVH #
        #######

        rvh_extractor = RoutingProductShapefileExtractor(
            routing_product_shp_path, hru_aspect_convention="ArcGIS"
        )
        rvh_config = rvh_extractor.extract()
        channel_profiles = rvh_config.pop("channel_profiles")

        gauge = [sb for sb in rvh_config["subbasins"] if sb.gauged]
        assert len(gauge) == 1
        gauge = gauge[0]

        rvh_config.update(
            {
                "land_subbasin_property_multiplier": SBGroupPropertyMultiplierCommand(
                    "Land", "MANNINGS_N", 1.0
                ),
                "lake_subbasin_property_multiplier": SBGroupPropertyMultiplierCommand(
                    "Lakes", "RESERVOIR_CREST_WIDTH", 1.0
                ),
            }
        )

        for k, v in rvh_config.items():
            model.config.rvh.update(k, v)

        #######
        # RVP #
        #######

        # The labels used for the following commands ("Lake_Soil_Lake_HRU", "Soil_Land_HRU", etc)
        # must correspond to the values of certain fields of the Routing Product:
        # LAND_USE_C, VEG_C, SOIL_PROF

        model.config.rvp.tmpl = """
        {soil_classes}

        {soil_profiles}

        {vegetation_classes}

        {land_use_classes}

        {avg_annual_runoff}

        {channel_profiles}
        """

        model.config.rvp.avg_annual_runoff = 594
        model.config.rvp.soil_classes = [SoilClassesCommand.Record("AQUIFER")]
        model.config.rvp.soil_profiles = [
            SoilProfilesCommand.Record("Lake_Soil_Lake_HRU", ("AQUIFER",), (5,)),
            SoilProfilesCommand.Record("Soil_Land_HRU", ("AQUIFER",), (5,)),
        ]
        model.config.rvp.vegetation_classes = [
            VegetationClassesCommand.Record("Veg_Land_HRU", 25, 5.0, 5.0),
            VegetationClassesCommand.Record("Veg_Lake_HRU", 0, 0, 0),
        ]
        model.config.rvp.land_use_classes = [
            LandUseClassesCommand.Record("Landuse_Land_HRU", 0, 1),
            LandUseClassesCommand.Record("Landuse_Lake_HRU", 0, 0),
        ]
        model.config.rvp.channel_profiles = channel_profiles

        #######
        # RVT #
        #######

        streaminputs_extractor = RoutingProductGridWeightExtractor(
            vic_streaminputs_nc_path, routing_product_shp_path
        )

        temperatures_extractor = RoutingProductGridWeightExtractor(
            vic_temperatures_nc_path, routing_product_shp_path
        )

        model.config.rvt.set_nc_variables(
            [
                dict(
                    name="StreamInputs",
                    data_type="PRECIP",
                    file_name_nc=vic_streaminputs_nc_path.name,
                    var_name_nc="Streaminputs",
                    dim_names_nc=("lon_dim", "lat_dim", "time"),
                    scale=4.0,
                    offset=0,
                    grid_weights=streaminputs_extractor.extract(),
                ),
                dict(
                    name="AverageTemp",
                    data_type="TEMP_AVE",
                    file_name_nc=vic_temperatures_nc_path.name,
                    var_name_nc="Avg_temp",
                    dim_names_nc=("lon_dim", "lat_dim", "time"),
                    grid_weights=temperatures_extractor.extract(),
                ),
                dict(
                    is_observation=True,
                    data_type="HYDROGRAPH",
                    subbasin_id=gauge.subbasin_id,
                    units="m3/s",
                    file_name_nc=observation_data_nc_path.name,
                    var_name_nc="Q",
                    dim_names_nc=("nstations", "time"),
                    index=1,  # StationIdx
                ),
            ]
        )

        #############
        # Run model #
        #############

        model(
            [
                vic_streaminputs_nc_path,
                vic_temperatures_nc_path,
                observation_data_nc_path,
            ],
        )

        ##########
        # Verify #
        ##########

        assert len(model.hydrograph.time) == (end - start).days + 1

        assert model.hydrograph.basin_name.item() == gauge.name

        csv_lines = model.outputs["diagnostics"].read_text().split("\n")
        assert csv_lines[1].split(",")[:-1] == [
            "HYDROGRAPH_ALL",
            observation_data_nc_path.name,
            "0.311049",
            "-11.9514",
            "0.471256",
        ]

        for d, q_sim in [
            (0, 85.97869699520872),
            (1000, 78.9516829670743),
            (2000, 70.29412760595676),
            (3000, 44.711237489482755),
            (4000, 129.98874279175033),
        ]:
            assert model.hydrograph.q_sim[d].item() == pytest.approx(q_sim)
