import datetime as dt

CF_RAVEN = {
    "TEMP_MIN": "tasmin",
    "TEMP_MAX": "tasmax",
    "TEMP_AVE": "tas",
    "RAINFALL": "rainfall",
    "PRECIP": "pr",
    "SNOWFALL": "prsn",
    "PET": "evspsbl",
    "HYDROGRAPH": "water_volume_transport_in_river_channel",
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
