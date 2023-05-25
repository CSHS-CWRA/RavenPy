import datetime as dt
import itertools
import re
from pathlib import Path
from textwrap import dedent, indent
from typing import (
    Any,
    Dict,
    Literal,
    Optional,
    Sequence,
    Tuple,
    Union,
    get_args,
    no_type_check,
)

import cftime  # noqa: F401
import xarray as xr
from pydantic import (
    Field,
    FilePath,
    HttpUrl,
    NonNegativeFloat,
    PositiveInt,
    PrivateAttr,
    ValidationError,
    confloat,
    validator,
)

from ..config import options
from .base import (
    Command,
    FlatCommand,
    GenericParameterList,
    ListCommand,
    ParameterList,
    Record,
    Sym,
)
from .utils import filter_for, nc_specs

"""
Notes
-----

TODO: Strip command outputs of white space to facilitate testing.
TODO: Add docstrings
TODO: Create tests in tests/config/test_commands.py

"""


INDENT = " " * 4
VALUE_PADDING = 10
T12 = Tuple[
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
    float,
]


# Validators
def reorder_time(cls, v):
    """TODO: Return dimensions as x, y, t. Currently only puts time at the end."""
    dims = list(v)
    for time_dim in ("t", "time"):
        if time_dim in dims:
            dims.remove(time_dim)
            dims.append(time_dim)
            break
    else:
        raise ValidationError("No time dimension found in dim_names_nc tuple.")

    return tuple(dims)


class Process(Command):
    """Process type embedded in HydrologicProcesses command.

    See processes.py for list of processes.
    """

    algo: str = "RAVEN_DEFAULT"
    source: str = None
    to: Sequence[str] = ()

    _sub = ""

    @validator("to", pre=True)
    def is_list(cls, v):
        if isinstance(v, str):
            return [
                v,
            ]
        return v

    @validator("source", "to", each_item=True)
    def is_state_variable(cls, v):
        if re.split(r"\[", v)[0] not in get_args(options.StateVariables):
            raise ValueError(f"{v} not recognized as a state variable.")
        return v

    def to_rv(self):
        t = " ".join(map("{:<19}".format, self.to))
        p = getattr(self, "p", None)
        p = p if p is not None else ""
        cmd = self._sub + self.__class__.__name__
        return f":{cmd:<20} {self.algo:<19} {self.source:<19} {t} {p}".strip()


class SoilModel(Command):
    __root__: int = None

    def to_rv(self):
        cmd = "SoilModel"
        return f":{cmd:<20} SOIL_MULTILAYER {self.__root__}\n"


class LinearTransform(Command):
    scale: float = 1
    offset: float = 0

    def to_rv(self):
        cmd = "LinearTransform"
        if (self.scale != 1) or (self.offset != 0):
            return f":{cmd:<20} {self.scale:.15f} {self.offset:.15f}\n"
        return ""


class RainSnowTransition(Command):
    """Specify the range of temperatures over which there will be a rain/snow mix when partitioning total precipitation into rain and snow components.

    Attributes
    ----------
    temp : float
      Midpoint of the temperature range [C].
    delta : float
      Range [C].
    """

    temp: Sym
    delta: Sym

    def to_rv(self):
        cmd = "RainSnowTransition"
        return f":{cmd:<20} {self.temp} {self.delta}\n"


class EvaluationPeriod(FlatCommand):
    """:EvaluationPeriod [period_name] [start yyyy-mm-dd] [end yyyy-mm-dd]"""

    name: str
    start: dt.date
    end: dt.date

    def to_rv(self):
        cmd = "EvaluationPeriod"
        return f":{cmd:<20} {self.name} {self.start} {self.end}\n"


class CustomOutput(FlatCommand):
    """
    Create custom output file to track a single variable, parameter or forcing function over time at a number of
    basins, HRUs, or across the watershed.

    Parameters
    ----------
    time_per : {'DAILY', 'MONTHLY', 'YEARLY', 'WATER_YEARLY', 'CONTINUOUS'}
      Time period.
    stat : {'AVERAGE', 'MAXIMUM', 'MINIMUM', 'RANGE', 'MEDIAN', 'QUARTILES', 'HISTOGRAM [min] [max] [# bins]'
      Statistic reported for each time inverval.
    variable: str
      Variable or parameter name. Consult the Raven documentation for the list of allowed names.
    space_agg : {'BY_BASIN', 'BY_HRU', 'BY_HRU_GROUP', 'BY_SB_GROUP', 'ENTIRE_WATERSHED'}
      Spatial evaluation domain.
    filename : str
      Output file name. Defaults to something approximately like `<run name>_<variable>_<time_per>_<stat>_<space_agg>.nc
    """

    time_per: Literal["DAILY", "MONTHLY", "YEARLY", "WATER_YEARLY", "CONTINUOUS"]
    stat: Literal["AVERAGE", "MAXIMUM", "MINIMUM", "RANGE", "MEDIAN", "QUARTILES"]
    variable: str
    space_agg: Literal[
        "BY_BASIN", "BY_HRU", "BY_HRU_GROUP", "BY_SB_GROUP", "ENTIRE_WATERSHED"
    ]
    filename: str = ""

    def to_rv(self):
        cmd = "CustomOutput"
        return f":{cmd:<20} {self.time_per} {self.stat} {self.variable} {self.space_agg} {self.filename}\n"


class SoilProfile(Record):
    name: str = ""
    soil_classes: Sequence[str] = ()
    thicknesses: Sequence[Sym] = ()

    def __str__(self):
        # From the Raven manual: {profile_name,#horizons,{soil_class_name,thick.}x{#horizons}}x[NP]
        n_horizons = len(self.soil_classes)
        horizon_data = list(itertools.chain(*zip(self.soil_classes, self.thicknesses)))
        fmt = "{:<16},{:>4}," + ",".join(n_horizons * ["{:>12}, {:>8}"])
        return fmt.format(self.name, n_horizons, *horizon_data)


class SoilProfiles(ListCommand):
    __root__: Sequence[SoilProfile]


class SubBasin(Record):
    """Record to populate RVH :SubBasins command internal table."""

    subbasin_id: PositiveInt = 1
    name: str = "sub_001"
    downstream_id: int = -1
    profile: str = "NONE"
    reach_length: float = 0
    gauged: bool = True
    gauge_id: Optional[str] = ""  # This attribute is not rendered to RVH

    def __str__(self):
        d = self.dict()
        d["reach_length"] = d["reach_length"]  # if d["reach_length"] else "ZERO-"
        d["gauged"] = int(d["gauged"])
        del d["gauge_id"]
        return " ".join(f"{v: <{VALUE_PADDING}}" for v in d.values())


class SubBasins(ListCommand):
    """SubBasins command (RVH)."""

    __root__: Sequence[SubBasin]

    _template = """
        :SubBasins
          :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
          :Units      none none          none    none           km    none
        {_records}
        :EndSubBasins
    """


class HRU(Record):
    """Record to populate :HRUs command internal table (RVH)."""

    hru_id: PositiveInt = 1
    area: Sym = 0  # km^2
    elevation: float = 0  # meters
    latitude: float = 0
    longitude: float = 0
    subbasin_id: PositiveInt = 1
    land_use_class: str = "[NONE]"
    veg_class: str = "[NONE]"
    soil_profile: str = "[NONE]"
    aquifer_profile: str = "[NONE]"
    terrain_class: str = "[NONE]"
    slope: NonNegativeFloat = 0.0
    aspect: confloat(ge=0, le=360) = 0.0
    # This field is not part of the Raven config, it is needed for serialization,
    # to specify which HRU subclass to use when necessary
    hru_type: Optional[str] = None

    def __str__(self):
        import numpy as np

        d = self.dict()
        del d["hru_type"]

        # Adjust horizontal spacing to match attributes names
        attrs = HRUs._template.splitlines()[2].strip()
        ns = np.array(list(map(len, re.split(r"\s[\w\:]", attrs))))
        ns[1:] += 1

        return " ".join(f"{v:<{n}}" for v, n in zip(d.values(), ns))


class HRUs(ListCommand):
    """HRUs command (RVH)."""

    __root__: Sequence[HRU]

    _template = """
        :HRUs
          :Attributes    AREA      ELEVATION       LATITUDE      LONGITUDE BASIN_ID       LAND_USE_CLASS         VEG_CLASS      SOIL_PROFILE  AQUIFER_PROFILE TERRAIN_CLASS      SLOPE     ASPECT
          :Units          km2              m            deg            deg     none                 none              none              none             none          none        deg       degN
        {_records}
        :EndHRUs
        """


class HRUGroup(FlatCommand):
    class _Rec(Record):
        __root__: Sequence[str]

        def __str__(self):
            return ",".join(self.__root__)

    name: str
    groups: _Rec

    _template = """
        :{_cmd} {name}
        {_records}
        :End{_cmd}
        """


class Reservoir(FlatCommand):
    """Reservoir command (RVH)."""

    name: str = "Lake_XXX"
    subbasin_id: int = Field(0, alias="SubBasinID")
    hru_id: int = Field(0, alias="HRUID")
    type: str = Field("RESROUTE_STANDARD", alias="Type")
    weir_coefficient: float = Field(0, alias="WeirCoefficient")
    crest_width: float = Field(0, alias="CrestWidth")
    max_depth: float = Field(0, alias="MaxDepth")
    lake_area: float = Field(0, alias="LakeArea", description="Lake area in m2")

    _template = """
            :Reservoir {name}
            {_commands}
            :EndReservoir
            """


class SubBasinGroup(FlatCommand):
    """SubBasinGroup command (RVH)."""

    name: str = ""
    sb_ids: Sequence[int] = ()

    def to_rv(self):
        template = """
            :SubBasinGroup {name}
              {sb_ids}
            :EndSubBasinGroup
            """
        d = self.dict()
        n_per_line = 10
        sbids = sorted(self.sb_ids)
        sbids_lines = [
            map(str, sbids[i : i + n_per_line])
            for i in range(0, len(sbids), n_per_line)
        ]
        d["sb_ids"] = "\n  ".join([", ".join(sbids) for sbids in sbids_lines])
        return dedent(template).format(**d)


class SBGroupPropertyMultiplier(FlatCommand):
    group_name: str
    parameter_name: str
    mult: float

    _template = "\n:SBGroupPropertyMultiplier {group_name} {parameter_name} {mult}"


# TODO: Convert to new config
class ChannelProfile(FlatCommand):
    """ChannelProfile command (RVP)."""

    name: str = "chn_XXX"
    bed_slope: float = 0
    survey_points: Tuple[Tuple[float, float], ...] = ()
    roughness_zones: Tuple[Tuple[float, float], ...] = ()

    def to_rv(self):
        template = """
            :ChannelProfile {name}
              :Bedslope {bed_slope}
              :SurveyPoints
            {survey_points}
              :EndSurveyPoints
              :RoughnessZones
            {roughness_zones}
              :EndRoughnessZones
            :EndChannelProfile
            """
        d = self.dict()
        d["survey_points"] = "\n".join(
            f"{INDENT * 2}{p[0]} {p[1]}" for p in d["survey_points"]
        )
        d["roughness_zones"] = "\n".join(
            f"{INDENT * 2}{z[0]} {z[1]}" for z in d["roughness_zones"]
        )
        return dedent(template).format(**d)


class GridWeights(Command):
    """GridWeights command.

    Important note: this command can be embedded in both a `GriddedForcing`
    or a `StationForcing`.

    The default is to have a single cell that covers an entire single HRU, with a
    weight of 1.
    """

    number_hrus: int = Field(1, alias="NumberHRUs")
    number_grid_cells: int = Field(1, alias="NumberGridCells")

    class GWRecord(Record):
        __root__: Tuple[int, int, float] = (1, 0, 1.0)

        def __iter__(self):
            return iter(self.__root__)

        def __getitem__(self, item):
            return self.__root__[item]

        def __str__(self):
            return " ".join(map(str, self.__root__))

    data: Sequence[GWRecord] = Field((GWRecord(),))

    @classmethod
    def parse(cls, s):
        pat = r"""
        :GridWeights
        \s*:NumberHRUs\s+(\d+)
        \s*:NumberGridCells\s+(\d+)
        \s*(.+)
        :EndGridWeights
        """
        m = re.match(dedent(pat).strip(), s.strip(), re.DOTALL)
        n_hrus, n_grid_cells, data = m.groups()  # type: ignore
        data = [d.strip().split() for d in data.split("\n")]
        data = tuple((int(h), int(c), float(w)) for h, c, w in data)
        return cls(
            number_hrus=int(n_hrus), number_grid_cells=int(n_grid_cells), data=data
        )


class RedirectToFile(Command):
    """RedirectToFile command (RVT).

    For the moment, this command can only be used in the context of a `GriddedForcingCommand`
    or a `StationForcingCommand`, as a `grid_weights` field replacement when inlining is not
    desired.
    """

    __root__: FilePath

    def to_rv(self):
        cmd = "RedirectToFile"
        return f":{cmd:<20} {self.__root__}\n"


class ReadFromNetCDF(FlatCommand):
    # TODO: When encoding Path, return only path.name
    # Order of HttpUrl, Path is important to avoid casting strings to Path.
    file_name_nc: Union[HttpUrl, Path] = Field(
        ..., alias="FileNameNC", description="NetCDF file name."
    )
    var_name_nc: str = Field(
        ..., alias="VarNameNC", description="NetCDF variable name."
    )
    dim_names_nc: Sequence[str] = Field(..., alias="DimNamesNC")
    station_idx: int = Field(
        1,
        alias="StationIdx",
        description="NetCDF index along station dimension. Starts at 1.",
    )
    time_shift: float = Field(
        None, alias="TimeShift", description="Time stamp shift in days."
    )
    linear_transform: LinearTransform = Field(None, alias="LinearTransform")
    deaccumulate: bool = Field(None, alias="Deaccumulate")

    latitude_var_name_nc: str = Field(None, alias="LatitudeVarNameNC")
    longitude_var_name_nc: str = Field(None, alias="LongitudeVarNameNC")
    elevation_var_name_nc: str = Field(None, alias="ElevationVarNameNC")

    _reorder_time = validator("dim_names_nc", allow_reuse=True)(reorder_time)

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=1, alt_names=(), **kwds):
        """Instantiate class from netCDF dataset."""
        specs = nc_specs(fn, data_type, station_idx, alt_names)
        specs.update(kwds)
        attrs = filter_for(cls, specs)
        return cls(**attrs)

    @property
    def da(self) -> xr.DataArray:
        """Return DataArray from configuration."""
        # TODO: Apply linear transform and time shift
        da = xr.open_dataset(self.file_name_nc)[self.var_name_nc]
        if len(self.dim_names_nc) == 1:
            return da
        elif len(self.dim_names_nc) == 2:
            return da.isel({self.dim_names_nc[0]: self.station_idx})
        else:
            raise NotImplementedError()


class GriddedForcing(ReadFromNetCDF):
    """GriddedForcing command (RVT)."""

    name: str = ""
    forcing_type: options.Forcings = Field(None, alias="ForcingType")
    grid_weights: Union[GridWeights, RedirectToFile] = Field(
        GridWeights(), alias="GridWeights"
    )
    # StationIdx is not relevant to GriddedForcing
    station_idx: int = None

    _template = """
    :{_cmd} {name}
    {_commands}
    :End{_cmd}
    """

    @validator("dim_names_nc")
    def check_dims(cls, v):
        if len(v) != 3:
            raise ValueError(
                "GriddedForcing netCDF datasets should have  three dimensions (lon, lat, time)."
            )
        return v


class StationForcing(GriddedForcing):
    """StationForcing command (RVT)."""

    @validator("dim_names_nc")
    def check_dims(cls, v):
        if len(v) != 2:
            raise ValueError(
                "StationForcing netCDF datasets should have two dimensions (station, time)."
            )
        return v

    @classmethod
    def from_nc(cls, fn, data_type, alt_names=(), **kwds):
        """Instantiate class from netCDF dataset."""
        specs = nc_specs(fn, data_type, station_idx=None, alt_names=alt_names)
        specs.update(kwds)

        attrs = filter_for(cls, specs)

        # Build default gridweights
        if "grid_weights" not in kwds and "GridWeights" not in kwds:
            size = specs["_dim_size_nc"]
            size.pop(specs["_time_dim_name_nc"])
            nstations = list(size.values())[0]
            data = [[i + 1, i + 1, 1] for i in range(nstations)]
            gw = GridWeights(
                number_hrus=nstations, number_grid_cells=nstations, data=data
            )
            attrs["grid_weights"] = gw

        return cls(**attrs)


class Data(FlatCommand):
    data_type: options.Forcings = ""
    units: str = ""
    read_from_netcdf: ReadFromNetCDF = Field(..., alias="ReadFromNetCDF")

    _template: str = """
                :{_cmd} {data_type} {units}
                {_commands}
                :End{_cmd}
                """

    _nested: bool = False

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=1, alt_names=(), **kwds):
        specs = nc_specs(fn, data_type, station_idx, alt_names)
        specs.update(kwds)
        attrs = filter_for(cls, specs)
        return cls(**attrs)


class Gauge(FlatCommand):
    name: str = "default"
    latitude: float = Field(..., alias="Latitude")
    longitude: float = Field(..., alias="Longitude")
    elevation: float = Field(None, alias="Elevation")

    rain_correction: Sym = Field(
        None, alias="RainCorrection", description="Rain correction"
    )
    snow_correction: Sym = Field(
        None, alias="SnowCorrection", description="Snow correction"
    )

    monthly_ave_evaporation: Sequence = Field(None, alias="MonthlyAveEvaporation")
    monthly_ave_temperature: Sequence = Field(None, alias="MonthlyAveTemperature")
    monthly_min_temperature: Sequence = Field(None, alias="MonthlyMinTemperature")
    monthly_max_temperature: Sequence = Field(None, alias="MonthlyMaxTemperature")

    data: Sequence[Data] = Field(None, alias="Data")

    _template = """
        :Gauge {name}
        {_commands}
        :EndGauge
        """

    _nested: bool = False

    @validator("monthly_ave_evaporation", "monthly_ave_temperature")
    def confirm_monthly(cls, v):
        if v is not None and len(v) != 12:
            raise ValidationError("One value per month needed.")
        return v

    @classmethod
    def from_nc(
        cls,
        fn: Union[Path, Sequence[Path]],
        data_type: Sequence[str] = None,
        station_idx: int = 1,
        alt_names: Optional[Dict[str, str]] = None,
        mon_ave: bool = False,
        data_kwds: Optional[
            Dict[
                str,
                Union[Dict[str, Union[bool, Dict[str, int], float]], Dict[str, Any]],
            ]
        ] = None,
        **kwds,
    ) -> "Gauge":
        """Return Gauge instance with configuration options inferred from the netCDF itself.

        Parameters
        ----------
        fn : Union[Path, Sequence[Path]],
          NetCDF file path or paths.
        data_type: Sequence[str], None
          Raven data types to extract from netCDF files, e.g. 'PRECIP', 'AVE_TEMP'. The algorithm tries to find all
          forcings in each file until one is found, then it stops searching for it in the following files.
        station_idx: int
          Index along station dimension. Starts at 1. Should be the same for all netCDF files.
        alt_names: dict
          Alternative variable names keyed by data type.
          Use this if variables do not correspond to CF standard defaults.
        mon_ave: bool
          If True, compute the monthly average.
        data_kwds: Dict[options.Forcings, Dict[str, str]]]
          Additional `:Data` parameters keyed by forcing type and station id. Overrides inferred parameters.
          Use keyword "ALL" to pass parameters to all variables.
        **kwds
          Additional arguments for Gauge.

        Returns
        -------
        Gauge
          Gauge instance.
        """
        if alt_names is None:
            alt_names = {}
        if data_kwds is None:
            data_kwds = {}

        if isinstance(fn, (str, Path)):
            fn = [fn]

        idx = station_idx
        data = {}
        attrs = {}

        forcings = set(data_type or get_args(options.Forcings))

        for f in fn:
            for dtype in forcings:
                try:
                    specs = nc_specs(
                        f, dtype, idx, alt_names.get(dtype, ()), mon_ave=mon_ave
                    )
                except ValueError:
                    pass

                else:
                    # Add extra keywords
                    ex_f = {**data_kwds.get(dtype, {}), **data_kwds.get("ALL", {})}

                    # Filter data attributes
                    data[dtype] = filter_for(Data, specs, **ex_f)

                    attrs.update(filter_for(cls, specs))

            # Remove forcings that have been found, so that they're not searched in the following files.
            forcings.difference_update(data)

        if len(data) == 0:
            raise ValueError("No data found in netCDF files.")

        # Default Gauge name
        attrs["name"] = attrs.get("name", f"Gauge_{idx}")

        attrs = filter_for(cls, attrs, **kwds, data=list(data.values()))
        return cls(**attrs)

    @property
    def ds(self) -> xr.Dataset:
        """Return xarray Dataset with forcing variables keyed by Raven forcing names."""
        ds = {}
        for data in self.data:
            ds[data.data_type] = data.read_from_netcdf.da
        return xr.Dataset(ds)


class ObservationData(Data):
    uid: str = "1"  # Basin ID or HRU ID
    data_type: Literal["HYDROGRAPH"] = "HYDROGRAPH"

    _template = """
        :ObservationData {data_type} {uid} {units}
        {_commands}
        :EndObservationData
        """

    @classmethod
    def from_nc(cls, fn, station_idx: int = 1, alt_names=(), **kwds):
        specs = nc_specs(fn, "HYDROGRAPH", station_idx, alt_names)
        attrs = filter_for(cls, specs, **kwds, data_type="HYDROGRAPH")
        return cls(**attrs)


class HRUState(Record):
    hru_id: int = 1
    data: Dict[str, Sym] = Field(default_factory=dict)

    def __str__(self):
        return ",".join(map(str, (self.hru_id,) + tuple(self.data.values())))

    @classmethod
    def parse(cls, sol, names=None):
        idx, *values = re.split(r",|\s+", sol.strip())
        if names is None:
            names = [f"S{i}" for i in range(len(values))]
        return cls(hru_id=idx, data=dict(zip(names, values)))


class HRUStateVariableTable(ListCommand):
    """Table of HRU state variables.

    If the HRUState include different attributes, the states will be modified to include all attributes.
    """

    __root__: Sequence[HRUState]

    _Attributes: Sequence[str] = PrivateAttr(None)

    def __init__(self, **data):
        from itertools import chain

        super().__init__(**data)

        # Get all attribute names from HRUStates
        names = (
            sorted(list(set(chain(*[tuple(s.data.keys()) for s in self.__root__]))))
            or None
        )

        self.__root__ = [
            HRUState(hru_id=s.hru_id, data={n: s.data.get(n, 0.0) for n in names})
            for s in self.__root__
        ]

        # Store names in private class attribute
        self._Attributes = names

    @classmethod
    def parse(cls, sol: str):
        from collections import defaultdict

        pat = re.compile(
            r"^:HRUStateVariableTable$\s*:Attributes,(?P<atts>.+)$\s*:Units,.*$\s*(?:("
            r"?P<vars>.+)$)+\s*:EndHRUStateVariableTable",
            re.MULTILINE,
        )
        matches = pat.findall(sol)

        d = defaultdict(dict)
        for m in matches:
            atts = [a.strip() for a in m[0].strip(",").split(",")]
            vals = [re.split(r",|\s+", line.strip(",")) for line in m[1:]]

            for val in vals:
                idx, *values = val
            d[idx].update(dict(zip(atts, values)))

        out = []
        for idx, data in d.items():
            out.append(HRUState(hru_id=idx, data=data))
        return cls(__root__=out)


class BasinIndex(Command):
    """Initial conditions for a flow segment."""

    # TODO: Check that qout and cie can be written with separating commas.
    sb_id: int = 1
    name: str = "watershed"
    channel_storage: float = Field(0.0, alias="ChannelStorage")
    rivulet_storage: float = Field(0.0, alias="RivuletStorage")
    qout: Sequence[float] = Field((1.0, 0.0, 0.0), alias="Qout")
    qlat: Sequence[float] = Field(None, alias="Qlat")
    qin: Sequence[float] = Field(None, alias="Qin")

    _template = """
         :BasinIndex {sb_id} {name}
         {_commands}
          """

    @classmethod
    @no_type_check
    def parse(cls, s):
        pat = r"""
        :BasinIndex (.+?)
        (.+)
        """
        m = re.search(dedent(pat).strip(), s, re.DOTALL)
        index_name = re.split(r",|\s+", m.group(1).strip())
        rec_values = {"sb_id": index_name[0], "name": index_name[1]}
        for line in m.group(2).strip().splitlines():
            all_values = filter(None, re.split(r",|\s+", line.strip()))
            cmd, *values = all_values
            if cmd in [":ChannelStorage", ":RivuletStorage"]:
                # FIXME: assertions should not be found outside of testing code. Replace with conditional logic.
                assert len(values) == 1
                rec_values[cmd[1:]] = float(values[0])
            else:
                rec_values[cmd[1:]] = tuple(values)
        return cls(**rec_values)


class BasinStateVariables(ListCommand):
    __root__: Sequence[BasinIndex]

    @classmethod
    @no_type_check
    def parse(cls, sol):
        pat = r"""
        :BasinStateVariables
        (.+)
        :EndBasinStateVariables
        """
        m = re.search(dedent(pat).strip(), sol, re.DOTALL)
        bis = filter(None, m.group(1).strip().split(":BasinIndex"))
        bs = [BasinIndex.parse(f":BasinIndex {bi}") for bi in bis]
        return cls.parse_obj(bs)


class SoilClasses(ListCommand):
    class SoilClass(Record):
        name: str
        mineral: Tuple[float, float, float] = None
        organic: float = None

        @validator("mineral")
        def validate_mineral(cls, v):
            """Assert sum of mineral fraction is 1."""
            if v is not None:
                # FIXME: assertions should not be found outside of testing code. Replace with conditional logic.
                assert sum(v) == 1, "Mineral fraction should sum to 1."
            return v

        @validator("mineral", "organic", each_item=True)
        def validate_pct(cls, v):
            if v is not None:
                # FIXME: assertions should not be found outside of testing code. Replace with conditional logic.
                assert (v >= 0) and (v <= 1), "Value should be in [0,1]."
            return v

        def __str__(self):
            if self.mineral is None and self.organic is None:
                return self.name
            else:
                fmt = "{name:<16},{mineral[0]:>18},{mineral[1]:>18},{mineral[2]:>18},{organic:>18}"
                return fmt.format(**self.dict())

    __root__: Sequence[SoilClass] = ()
    _Attributes: Sequence[str] = PrivateAttr(["%SAND", "%CLAY", "%SILT", "%ORGANIC"])
    _Units: Sequence[str] = PrivateAttr(["none", "none", "none", "none"])


class VegetationClass(Record):
    name: str = ""
    max_ht: float = 0.0
    max_lai: float = 0.0
    max_leaf_cond: float = 0.0

    def __str__(self):
        template = "{name:<16},{max_ht:>18},{max_lai:>18},{max_leaf_cond:>18}"
        return template.format(**self.dict())


class VegetationClasses(ListCommand):
    __root__: Sequence[VegetationClass]
    _Attributes: Sequence[str] = PrivateAttr(["MAX_HT", "MAX_LAI", "MAX_LEAF_COND"])
    _Units: Sequence[str] = PrivateAttr(["m", "none", "mm_per_s"])


class LandUseClass(Record):
    name: str = ""
    impermeable_frac: Sym = 0.0
    forest_coverage: Sym = 0.0

    def __str__(self):
        template = "{name:<16},{impermeable_frac:>18},{forest_coverage:>18}"
        return template.format(**self.dict())


class LandUseClasses(ListCommand):
    __root__: Sequence[LandUseClass]
    _Attributes: Sequence[str] = PrivateAttr(["IMPERMEABLE_FRAC", "FOREST_COVERAGE"])
    _Units: Sequence[str] = PrivateAttr(["fract", "fract"])


class TerrainClass(Record):
    name: str
    hillslope_length: Sym
    drainage_density: Sym
    topmodel_lambda: Sym = None

    def __str__(self):
        out = f"{self.name:<16},{self.hillslope_length:>18},{self.drainage_density:>18}"
        if self.topmodel_lambda is not None:
            out += f",{self.topmodel_lambda:>18}"
        return out


class TerrainClasses(ListCommand):
    __root__: Sequence[TerrainClass]
    _Attributes: Sequence[str] = PrivateAttr(
        ["HILLSLOPE_LENGTH", "DRAINAGE_DENSITY", "TOPMODEL_LAMBDA"]
    )
    _Units: Sequence[str] = PrivateAttr(["m", "km/km2"])


class _MonthlyRecord:
    name: str = "[DEFAULT]"
    values: Sequence[float] = 12 * [
        1.0,
    ]

    def __str__(self):
        return f"{self.name:<16} " + ", ".join(map(str, self.values))


class SeasonalRelativeLAI(ListCommand):
    __root__: Sequence[_MonthlyRecord] = (_MonthlyRecord(),)


class SeasonalRelativeHeight(ListCommand):
    __root__: Sequence[_MonthlyRecord] = (_MonthlyRecord(),)


class SubBasinProperty(Record):
    sb_id: str
    values: Sequence[Sym]

    def __str__(self):
        return indent(
            ", ".join(
                [
                    self.sb_id,
                ]
                + list(map(str, self.values))
            ),
            Command._indent,
        )


class SubBasinProperties(Command):
    parameters: Sequence[options.SubBasinProperties] = Field(None, alias="Parameters")
    records: Sequence[SubBasinProperty] = Field(None)


class SoilParameterList(GenericParameterList):
    parameters: Sequence[options.SoilParameters] = Field(None, alias="Parameters")


class VegetationParameterList(GenericParameterList):
    parameters: Sequence[options.VegetationParameters] = Field(None, alias="Parameters")


class LandUseParameterList(GenericParameterList):
    parameters: Sequence[options.LandUseParameters] = Field(None, alias="Parameters")


class EnsembleMode(Command):
    n: int = Field(..., description="Number of members")
    mode: Literal["ENSEMBLE_ENKF"] = "ENSEMBLE_ENKF"

    _template = ":{_cmd} {mode} {n}\n"


class ObservationErrorModel(FlatCommand):
    state: Literal["STREAMFLOW"]
    dist: Literal["DIST_UNIFORM", "DIST_NORMAL", "DIST_GAMMA"]
    p1: float
    p2: float
    adj: Literal["ADDITIVE", "MULTIPLICATIVE"]

    _template = ":{_cmd} {state} {dist} {p1} {p2} {adj}\n"


class ForcingPerturbation(FlatCommand):
    forcing: options.Forcings
    dist: Literal["DIST_UNIFORM", "DIST_NORMAL", "DIST_GAMMA"]
    p1: float
    p2: float
    adj: Literal["ADDITIVE", "MULTIPLICATIVE"]
    hru_grp: str = ""

    _template = ":{_cmd} {forcing} {dist} {p1} {p2} {adj} {hru_grp}\n"


class AssimilatedState(FlatCommand):
    state: Union[options.StateVariables, Literal["STREAMFLOW"]]
    group: str

    _template = ":{_cmd} {state} {group}\n"


class AssimilateStreamflow(FlatCommand):
    sb_id: str

    _template = ":{_cmd} {sb_id}\n"


# TODO: Harmonize this...
# Aliases for convenience
SC = SoilClasses.SoilClass
LU = LandUseClass
VC = VegetationClass
TC = TerrainClass
SB = SubBasin
SP = SoilProfile
PL = ParameterList


# ---
