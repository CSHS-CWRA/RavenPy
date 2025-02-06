# Taken from https://cfconventions.org/Data/cf-standard-names/current/build/cf-standard-name-table.html
# Augmented with https://pcmdi.llnl.gov/mips/cmip3/variableList.html
CF_RAVEN = {
    "TEMP_MIN": ("tasmin",),
    "TEMP_DAILY_MIN": ("tasmin",),
    "TEMP_MAX": ("tasmax",),
    "TEMP_DAILY_MAX": ("tasmax",),
    "TEMP_AVE": ("tas", "air_temperature"),
    "TEMP_DAILY_AVE": ("tas", "air_temperature"),
    "RAINFALL": (
        "rainfall_flux",
        "prra",
        "rainfall",
    ),
    "PRECIP": ("pr", "precipitation_flux"),
    "PRECIP_DAILY_AVE": ("pr", "precipitation_flux"),
    "SNOWFALL": ("prsn", "snowfall_flux", "snowfall"),
    "REL_HUMIDITY": ("hurs", "relative_humidity"),
    "AIR_PRES": ("ps", "surface_air_pressure"),
    "SHORTWAVE": ("rsds", "surface_downwelling_shortwave_flux_in_air"),
    "LW_INCOMING": ("rlds", "surface_downwelling_longwave_flux_in_air"),
    "CLOUD_COVER": ("clt", "cloud_area_fraction"),
    "WIND_VEL": (
        "sfcWind",
        "wind_speed",
    ),
    "PET": ("evspsbl", "water_evapotranspiration_flux", "water_evaporation_flux"),
    "HYDROGRAPH": ("water_volume_transport_in_river_channel",),
}


NetCDFAttribute = {
    "title": "Simulated river discharge",
    "history": "Created on {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} by Raven",
    "references": "Craig, J.R., and the Raven Development Team, Raven user's and developer's manual ("
    "Version 2.8), URL: http://raven.uwaterloo.ca/ (2018).",
    "comment": "Raven Hydrological Framework version {raven_version}",
    "model_id": "{identifier}",
    "time_frequency": "day",
    "time_coverage_start": "{start_date}",
    "time_coverage_end": "{end_date}",
}

MonthlyAverages = {
    "TEMP_AVE": "MonthlyAveTemperature",
    "TEMP_MIN": "MonthlyMinTemperature",
    "TEMP_MAX": "MonthlyMaxTemperature",
    "PET": "MonthlyAveEvaporation",
}

RAVEN_OUTPUT_FMT = {
    "solution": "{run_name}solution.rvc",
    "hydrograph": "{run_name}Hydrographs.nc",
    "storage": "{run_name}WatershedStorage.nc",
    "diagnostics": "{run_name}Diagnostics.csv",
    "messages": "Raven_errors.txt",
}
