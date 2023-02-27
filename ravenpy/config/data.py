from pydantic.dataclasses import dataclass
from pydantic import validator, parse_obj_as
from base import RavenCommand, RavenValue, RavenSwitch, Command, RavenList, RavenOption, CF_RAVEN
from typing import Union, List, Tuple
from pathlib import Path
from dataclasses import asdict
import xarray as xr
from textwrap import dedent
import options


@dataclass
class FileNameNc(RavenValue):
    """NetCDF file name."""
    # Can be an url or a path; note that the order "str, Path" is important here because
    # otherwise pydantic would try to coerce a string into a Path, which is a problem for a url
    # because it messes with slashes; so a Path will be one only when explicitly specified.
    value: Union[str, Path] = ""

@dataclass
class VarNameNc(RavenValue):
    """NetCDF variable name."""


@dataclass
class DimNamesNC(RavenList):
    value: List[str]

    @validator("value", pre=True)
    def time_is_last(cls, v):
        """Return dimensions with time as the last dimension."""
        v = list(v)
        for time_dim in ("t", "time"):
            if time_dim in v:
                v.remove(time_dim)
                v.append(time_dim)
                break
        else:
            raise ValueError("No time dimension found.")
        return v

@dataclass
class StationIdx(RavenValue):
    """NetCDF index along station dimension. Starts at 1."""
    value: Union[int, str] = 1


@dataclass
class TimeShift(RavenValue):
    """Time stamp shift in days."""
    value: float


@dataclass
class LinearTransform(RavenCommand):
    scale: float = 1
    offset: float = 0

    def to_rv(self):
        template = ":LinearTransform {scale:.15f} {offset:.15f}\n"
        if (self.scale != 1) or (self.offset != 0):
            return template.format(**asdict(self))
        return ""

@dataclass
class Deaccumulate(RavenSwitch):
    """Deaccumulate variable."""


@dataclass
class Latitude(RavenCommand):
    value: float


@dataclass
class Longitude(RavenCommand):
    value: float


@dataclass
class Elevation(RavenCommand):
    value: float

@dataclass
class RainCorrection(RavenValue):
    """Rain correction"""
    value: float

@dataclass
class SnowCorrection(RavenValue):
    """Snow correction"""
    value: float


@dataclass
class MonthlyAveEvaporation(RavenList):
    value: List[float]

    @validator('value')
    def has_12_months(cls, v):
        if len(v) != 12:
            raise ValueError("There should be 12 monthly values.")


@dataclass
class MonthlyAveTemperature(MonthlyAveEvaporation):
    """Monthly average temperature."""

@dataclass
class ForcingType(RavenOption):
    value: options.Forcings


@dataclass
class LatitudeVarNameNC(RavenValue):
    """"""

@dataclass
class LongitudeVarNameNC(RavenValue):
    """"""

@dataclass
class ElevationVarNameNC(RavenValue):
    """"""

@dataclass
class NumberHRUs(RavenValue):
    value: int

@dataclass
class NumberGridCells(RavenValue):
    value: int


@dataclass
class GridWeights(RavenCommand):
    """GridWeights command.

    Important note: this command can be embedded in both a `GriddedForcingCommand` or a `StationForcingCommand`.
    The default is to have a single cell that covers an entire single HRU, with a weight of 1.
    """

    number_hrus: NumberHRUs = 1
    number_grid_cells: NumberGridCells = 1
    data: List[Tuple[int, int, float]] = ((1, 0, 1.0),)

    @classmethod
    def parse(cls, s):
        import re
        pat = r"""
        :GridWeights
          :NumberHRUs (\d+)
          :NumberGridCells (\d+)
          (.+)
        :EndGridWeights
        """
        m = re.match(dedent(pat).strip(), s, re.DOTALL)
        n_hrus, n_grid_cells, data = m.groups()  # type: ignore
        data = [d.strip().split() for d in data.split("\n")]
        data = tuple((int(h), int(c), float(w)) for h, c, w in data)
        return cls(
            number_hrus=int(n_hrus), number_grid_cells=int(n_grid_cells), data=data
        )

    def to_rv(self):
        template = """
        :GridWeights
          {number_hrus}
          {number_grid_cells}
        {data}
        :EndGridWeights
        """

        d = asdict(self)
        d["data"] = "\n".join(f"  {p[0]} {p[1]} {p[2]}" for p in self.data)
        return dedent(template).strip().format(**d)



@dataclass
class ReadFromNetCDF(Command):
    file_name_nc: FileNameNc
    var_name_nc: VarNameNc
    dim_names_nc: DimNamesNC
    station_idx: StationIdx = StationIdx()
    time_shift: TimeShift = None
    linear_transform: LinearTransform = None
    deaccumulate: Deaccumulate = None

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=1, alt_names=()):
        """Instantiate class from netCDF dataset."""
        specs = nc_specs(fn, data_type, station_idx, alt_names)
        attrs = filter_for(cls, specs)
        return cls(station_idx=StationIdx(station_idx),
                   **attrs)


@dataclass
class Data(Command):
    data_type: options.Forcings = ""
    units: str = ""
    read_from_netcdf: ReadFromNetCDF = None

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=1, alt_names=()):
        specs = nc_specs(fn, data_type, station_idx, alt_names)
        return cls(data_type=data_type,
                   units=specs.pop("units", None),
                   read_from_netcdf=filter_for(ReadFromNetCDF, specs))

@dataclass
class ObservationData(Data):
    data_type: str = None
    id: str = None  # HRU_ID or SUBBASIN_ID
    units: str = None
    read_from_netcdf: ReadFromNetCDF = None

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=1, alt_names=()):
        specs = nc_specs(fn, data_type, station_idx, alt_names)
        return cls(data_type=data_type,
                   id=id,
                   units=specs.pop("units", None),
                   read_from_netcdf=filter_for(ReadFromNetCDF, specs))


@dataclass
class Gauge(Command):
    """One gauge includes multiple Data commands."""
    name: str
    latitude: Latitude
    longitude: Longitude
    elevation: Elevation
    rain_correction: RainCorrection = None
    snow_correction: SnowCorrection = None
    monthly_ave_evaporation: MonthlyAveEvaporation = None
    monthly_ave_temperature: MonthlyAveTemperature = None
    data: List[Data] = None

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=(1,), alt_names=()):
        for idx in station_idx:
            specs = nc_specs(fn, data_type, station_idx, alt_names)


@dataclass
class StationForcing(ReadFromNetCDF):
    forcing_name: str = None
    forcing_type: ForcingType = None
    grid_weights: GridWeights = None
    #map_station_to: MapStationTo
    latitude_var_name_nc: LatitudeVarNameNC = None
    longitude_var_name_nc: LongitudeVarNameNC = None
    elevation_var_name_nc: ElevationVarNameNC = None

    @validator("dim_names_nc")
    def check_2_dims(cls, v):
        if len(v) != 2:
            raise ValueError("StationForcing netCDF datasets should have two dimensions (station, time).")


@dataclass
class GriddedForcing(StationForcing):
    """GriddedForcing command (RVT)."""

    @validator("dim_names_nc")
    def check_3_dims(cls, v):
        if len(v) != 3:
            raise ValueError("GriddedForcing netCDF datasets should have three dimensions (lon, lat, time).")



def nc_specs(fn, data_type, station_idx, alt_names=()):
    """Extract specifications from netCDF file.

    Parameters
    ----------
    fn : str, Path
      NetCDF file path.
    data_type: str
      Raven data type.
    station_idx: int
      Index along station dimension. Starts at 1.
    alt_names: list
      Alternative variable names for data type if not the CF standard default.
    """
    from ravenpy.utilities.coords import infer_scale_and_offset

    # Convert to NumPy 0-based indexing
    i = station_idx - 1

    attrs = {"file_name_nc": {"value": fn},
             "data_type": data_type}
    with xr.open_dataset(fn) as ds:
        var_names = [CF_RAVEN[data_type], ] + list(alt_names)
        for v in var_names:
            if v in ds.data_vars:
                nc_var = ds[v]
                attrs["var_name_nc"] = {"value": v}
                attrs["dim_names_nc"] = {"value": nc_var.dims}
                attrs["units"] = nc_var.attrs.get("units")
                if attrs["units"] is not None:
                    s, o = infer_scale_and_offset(nc_var, data_type)
                    attrs["linear_transform"] = dict(scale=s, offset=o)

                break
        else:
            raise ValueError(f"No variable found for {data_type}.\n {ds.data_vars}")

        try:
            attrs["latitude_var_name_nc"] = {"value": ds.cf["latitude"].name}
            attrs["longitude_var_name_nc"] = {"value": ds.cf["longitude"].name}

            attrs["latitude"] = {"value": ds.cf["latitude"][i]}
            attrs["longitude"] = {"value": ds.cf["longitude"][i]}

        except KeyError:
            pass

        try:
            nc_elev = ds.cf["vertical"].name
            attrs["elevation"] = {"value": ds.cf["vertical"][i]}
        except KeyError:
            nc_elev = "elevation" if "elevation" in ds else None
        finally:
            if nc_elev is not None:
                attrs["elevation_var_name_nc"] = {"value": nc_elev}
                attrs["elevation"] = {"value": ds["elevation"][i]}

        if "station_id" in ds:
            if ds["station_id"].shape and len(ds["station_id"]) > i:
                attrs["name"] = {"value": ds["station_id"].values[i]}

        return attrs


def filter_for(kls, attrs):
    """Return attributes that are fields of dataclass."""
    return {k: v for (k,v) in attrs.items() if k in kls.__dataclass_fields__.keys()}
