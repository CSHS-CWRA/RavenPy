from .regionalization import read_gauged_params, read_gauged_properties, regionalize

dev_import_error_message = (
    "`{}` requires installation of the RavenPy development libraries. These can be installed using the"
    " `pip install ravenpy[dev]` recipe or via Anaconda (`conda env update -n ravenpy-env -f environment.yml`)"
    " from the RavenPy repository source files."
)
