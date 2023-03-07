def test_nc_specs(get_local_testdata):
    from ravenpy.new_config.utils import nc_specs

    f = get_local_testdata(
        "raven-gr4j-cemaneige/Salmon-River-Near-Prince-George_meteo_daily.nc"
    )
    attrs = nc_specs(f, "PRECIP", station_idx=1, alt_names=("rain",))
    assert "file_name_nc" in attrs
