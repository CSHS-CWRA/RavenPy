units = {
    "PRECIP": "mm/d",
    "PRECIP_DAILY_AVE": "mm/d",
    "PRECIP_5DAY": "mm",
    "SNOW_FRAC": "",
    "SNOWFALL": "mm/d",
    "RAINFALL": "mm/d",
    "RECHARGE": "mm/d",
    "TEMP_AVE": "degC",
    "TEMP_DAILY_AVE": "degC",
    "TEMP_MIN": "degC",
    "TEMP_DAILY_MIN": "degC",
    "TEMP_MAX": "degC",
    "TEMP_DAILY_MAX": "degC",
    "TEMP_MONTH_MAX": "degC",
    "TEMP_MONTH_MIN": "degC",
    "TEMP_MONTH_AVE": "degC",
    "TEMP_AVE_UNC": "degC",
    "TEMP_MAX_UNC": "degC",
    "TEMP_MIN_UNC": "degC",
    "AIR_DENS": "kg/m**3",
    "AIR_PRES": "kPa",
    "REL_HUMIDITY": "",
    "ET_RADIA": "MJ/m**2/d",
    "SHORTWAVE": "MJ/m**2/d",
    "SW_RADIA": "MJ/m**2/d",
    "SW_RADIA_NET": "MJ/m**2/d",
    "LW_RADIA_NET": "MJ/m**2/d",
    "LW_INCOMING": "MJ/m**2/d",
    "CLOUD_COVER": "",
    "DAY_LENGTH": "d",
    "DAY_ANGLE": "",
    "WIND_VEL": "m/s",
    "PET": "mm/d",
    "OW_PET": "mm/d",
    "PET_MONTH_AVE": "mm/d",
    "POTENTIAL_MELT": "mm/d",
    "SUBDAILY_CORR": "",
    "HYDROGRAPH": "m**3/s",
}

RAVEN_NO_DATA_VALUE = -1.2345
CALENDAR = "PROLEPTIC_GREGORIAN"


def nc_attrs(cls, val):
    """Ensure default netCDF attributes are present."""
    # FIXME: assertions should not be found outside of testing code. Replace with conditional logic.
    assert "model_id" in val  # noqa: S101

    out = default_nc_attrs()
    out.update(val)
    return out


def default_nc_attrs():
    """Return default NetCDF global attributes."""
    import datetime as dt

    now = dt.datetime.now().isoformat(timespec="seconds")
    # TODO: get version from shared library
    version = "3.7"

    return {
        "history": f"Created on {now} by Raven {version}",
        "references": "Craig, J.R., and the Raven Development Team, Raven user's and developer's manual "
        f"(Version {version}), URL: http://raven.uwaterloo.ca/ (2023).",
    }
