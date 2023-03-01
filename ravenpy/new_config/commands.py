import datetime as dt
import itertools
import re
from dataclasses import asdict, field
from itertools import chain
from pathlib import Path
from textwrap import dedent, indent
from typing import (
    Any,
    ClassVar,
    Dict,
    List,
    Literal,
    Optional,
    Sequence,
    Tuple,
    Union,
    no_type_check,
)

from pydantic import (
    BaseModel,
    Field,
    HttpUrl,
    PrivateAttr,
    ValidationError,
    root_validator,
    validator,
)
from pydantic.dataclasses import dataclass

from ravenpy.config import options

from .base import Command, ParameterList, ParameterListCommand, Sym
from .data import filter_for, nc_specs

INDENT = " " * 4
VALUE_PADDING = 10


# Validators
def reorder_time(cls, v):
    """Return dimensions with time as the last dimension."""
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
    source: list = ()
    to: list = ()

    def to_rv(self):
        s = " ".join(self.source)
        t = ("To" + " ".join(self.to)) if self.to else " "
        return f":{self.__class__.__name__} {s} {t}"


class SoilModel(Command):
    n: int = None

    def to_rv(self):
        cmd = "SoilModel"
        return f":{cmd:<20} SOIL_MULTILAYER {self.n}\n"


class LinearTransform(Command):
    scale: float = 1
    offset: float = 0

    def to_rv(self):
        if (self.scale != 1) or (self.offset != 0):
            return f":LinearTransform {self.scale:.15f} {self.offset:.15f}\n"
        return ""


class RainSnowTransition(Command):
    """Specify the range of temperatures over which there will be a rain/snow mix when partitioning total
    precipitation into rain and snow components.

    Attributes
    ----------
    temp : float
      Midpoint of the temperature range [C].
    delta : float
      Range [C].
    """

    temp: float
    delta: float

    def to_rv(self):
        return f":RainSnowTransition {self.temp} {self.delta}\n"


class EvaluationPeriod(Command):
    """:EvaluationPeriod [period_name] [start yyyy-mm-dd] [end yyyy-mm-dd]"""

    name: str
    start: dt.date
    end: dt.date

    def to_rv(self):
        return f":EvaluationPeriod {self.name} {self.start} {self.end}"


class CustomOutput(Command):
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
        template = ":CustomOutput {time_per} {stat} {variable} {space_agg} {filename}\n"
        return template.format(**self.dict())


class SoilProfile(Command):
    name: str = ""
    soil_classes: Sequence[str] = ()
    thicknesses: Sequence[Sym] = ()

    def to_rv(self):
        # From the Raven manual: {profile_name,#horizons,{soil_class_name,thick.}x{#horizons}}x[NP]
        n_horizons = len(self.soil_classes)
        horizon_data = list(itertools.chain(*zip(self.soil_classes, self.thicknesses)))
        fmt = "{:<16},{:>4}," + ",".join(n_horizons * ["{:>12},{:>6}"])
        return fmt.format(self.name, n_horizons, *horizon_data)


class SubBasin(Command):
    """Record to populate RVH :SubBasins command internal table."""

    subbasin_id: int = 1
    name: str = "sub_001"
    downstream_id: int = -1
    profile: str = "None"
    reach_length: float = 0
    gauged: bool = True
    gauge_id: Optional[str] = ""  # This attribute is not rendered to RVH

    def to_rv(self):
        d = self.dict()
        d["reach_length"] = d["reach_length"] if d["reach_length"] else "ZERO-"
        d["gauged"] = int(d["gauged"])
        del d["gauge_id"]
        return " ".join(f"{v: <{VALUE_PADDING}}" for v in d.values())


class SubBasins(Command):
    """SubBasins command (RVH)."""

    __root__: Sequence[SubBasin]

    _template = """
        :SubBasins
          :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
          :Units      none none          none    none           km    none
        {_commands}
        :EndSubBasins
    """


class HRU(Command):
    """Record to populate :HRUs command internal table (RVH)."""

    hru_id: int = 1
    area: float = 0  # km^2
    elevation: float = 0  # meters
    latitude: float = 0
    longitude: float = 0
    subbasin_id: int = 1
    land_use_class: str = ""
    veg_class: str = ""
    soil_profile: str = ""
    aquifer_profile: str = ""
    terrain_class: str = ""
    slope: float = 0.0
    aspect: float = 0.0
    # This field is not part of the Raven config, it is needed for serialization,
    # to specify which HRU subclass to use when necessary
    hru_type: Optional[str] = None

    def to_rv(self):
        d = self.dict()
        del d["hru_type"]
        return " ".join(f"{v: <{VALUE_PADDING * 2}}" for v in d.values())


class HRUsCommand(Command):
    """HRUs command (RVH)."""

    __root__: Sequence[HRU]

    _template = """
        :HRUs
          :Attributes      AREA  ELEVATION       LATITUDE      LONGITUDE BASIN_ID       LAND_USE_CLASS            VEG_CLASS      SOIL_PROFILE  AQUIFER_PROFILE TERRAIN_CLASS      SLOPE     ASPECT
          :Units            km2          m            deg            deg     none                  none                none              none             none          none        deg       degN
        {_commands}
        :EndHRUs
        """


class Reservoir(Command):
    """Reservoir command (RVH)."""

    subbasin_id: int = 0
    hru_id: int = 0
    name: str = "Lake_XXX"
    weir_coefficient: float = 0
    crest_width: float = 0
    max_depth: float = 0
    lake_area: float = 0  # in m^2

    def to_rv(self):
        template = """
            :Reservoir {name}
              :SubBasinID {subbasin_id}
              :HRUID {hru_id}
              :Type RESROUTE_STANDARD
              :WeirCoefficient {weir_coefficient}
              :CrestWidth {crest_width}
              :MaxDepth {max_depth}
              :LakeArea {lake_area}
            :EndReservoir
            """
        return dedent(template).format(**self.dict())


class SubBasinGroup(Command):
    """SubBasinGroup command (RVH)."""

    name: str = ""
    sb_ids: Tuple[int, ...] = ()

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


class SBGroupPropertyMultiplierCommand(Command):
    group_name: str
    parameter_name: str
    mult: float

    def to_rv(self):
        template = ":SBGroupPropertyMultiplier {group_name} {parameter_name} {mult}"
        return dedent(template).format(**self.dict())


class ChannelProfile(Command):
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
    data: Tuple[Tuple[int, int, float], ...] = ((1, 0, 1.0),)

    data_str: str = None  # Will be created in the __init__

    _template = """
            :GridWeights
            {_commands}
            {data_str}
            :EndGridWeights
            """

    def __init__(self, **data):
        super().__init__(**data)
        self.data_str = indent(
            "\n".join(f"{p[0]} {p[1]} {p[2]}" for p in self.data), "  "
        )

    @classmethod
    def parse(cls, s):
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


class RedirectToFile(Command):
    """RedirectToFile command (RVT).

    For the moment, this command can only be used in the context of a `GriddedForcingCommand`
    or a `StationForcingCommand`, as a `grid_weights` field replacement when inlining is not
    desired.
    """

    path: Path

    _template = ":{_cmd} {path}"

    @validator("path")
    def check_exists(cls, v):
        if not v.exists():
            raise ValidationError(f"File not found {v}")
        return v


class ReadFromNetCDF(Command):
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
    longitude_var_name_nc: str = Field(None, alias="LonitudeVarNameNC")
    elevation_var_name_nc: str = Field(None, alias="ElevationVarNameNC")

    _reorder_time = validator("dim_names_nc", allow_reuse=True)(reorder_time)

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=1, alt_names=()):
        """Instantiate class from netCDF dataset."""
        specs = nc_specs(fn, data_type, station_idx, alt_names)
        attrs = filter_for(cls, specs)
        return cls(station_idx=station_idx, **attrs)


class GriddedForcing(ReadFromNetCDF):
    """GriddedForcing command (RVT)."""

    name: str = ""
    forcing_type: options.Forcings = Field(None, alias="ForcingType")
    grid_weights: Union[GridWeights, RedirectToFile] = Field(
        GridWeights(), alias="GridWeights"
    )

    _template = """
    :{_cmd} {name}
    {_commands}
    :End{_cmd}
    """

    @validator("dim_names_nc")
    def check_3_dims(cls, v):
        if len(v) != 3:
            raise ValueError(
                "GriddedForcing netCDF datasets should have three dimensions (lon, lat, time)."
            )


class StationForcing(GriddedForcing):
    """StationForcing command (RVT)."""

    @validator("dim_names_nc")
    def check_2_dims(cls, v):
        if len(v) != 2:
            raise ValueError(
                "StationForcing netCDF datasets should have two dimensions (station, time)."
            )


class Data(Command):
    data_type: options.Forcings = ""
    units: str = ""
    read_from_netcdf: ReadFromNetCDF = Field(..., alias="ReadFromNetCDF")

    _template: str = """
                :{_cmd} {data_type} {units}
                {_commands}
                :End{_cmd}
                """

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=1, alt_names=()):
        specs = nc_specs(fn, data_type, station_idx, alt_names)
        return cls(
            data_type=data_type,
            units=specs.pop("units", None),
            read_from_netcdf=filter_for(ReadFromNetCDF, specs),
        )


class Gauge(Command):
    name: str = "default"
    latitude: float = Field(0, alias="Latitude")
    longitude: float = Field(0, alias="Lonitude")
    elevation: float = Field(0, alias="Elevation")

    # Accept strings to embed parameter names into Ostrich templates
    rain_correction: float = Field(
        None, alias="RainCorrection", description="Rain correction"
    )
    snow_correction: float = Field(
        None, alias="SnowCorrection", description="Snow correction"
    )

    monthly_ave_evaporation: Tuple[float] = Field(None, alias="MonthlyAveEvaporation")
    monthly_ave_temperature: Tuple[float] = Field(None, alias="MonthlyAveTemperature")

    data: Sequence[Data] = Field(None, alias="Data")

    _template = """
        :Gauge {name}
        {_commands}
        :EndGauge
        """

    @validator("monthly_ave_evaporation", "monthly_ave_temperature")
    def confirm_monthly(cls, v):
        if len(v) != 12:
            raise ValidationError("One value per month needed.")
        return v

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=(1,), alt_names=()):
        for idx in station_idx:
            specs = nc_specs(fn, data_type, station_idx, alt_names)



class ObservationData(Data):
    uid: int = 1  # Basin ID or HRU ID

    _template = """
        :ObservationData {data_type} {uid} {units}
        {_commands}
        :EndObservationData
        """

    @classmethod
    def from_nc(cls, fn, data_type, station_idx=1, uid=1, alt_names=()):
        specs = nc_specs(fn, data_type, station_idx, alt_names)
        return cls(
            data_type=data_type,
            uid=uid,
            units=specs.pop("units", None),
            read_from_netcdf=filter_for(ReadFromNetCDF, specs),
        )


class HRUState(Command):
    index: int = 1
    data: Dict[str, float] = Field(default_factory=dict)

    def to_rv(self):
        return ",".join(map(str, (self.index,) + tuple(self.data.values())))

    @classmethod
    def parse(cls, sol, names=None):
        idx, *values = re.split(r",|\s+", sol.strip())
        if names is None:
            names = [f"S{i}" for i in range(len(values))]
        return cls(index=idx, data=dict(zip(names, values)))


class HRUStateVariableTable(Command):
    __root__: List[HRUState]

    _attributes_str: Sequence[str] = PrivateAttr(None)

    _template = """
                 :HRUStateVariableTable
                   :Attributes,{_attributes_str}
                   :Units
                 {_commands}
                 :EndHRUStateVariableTable
                 """

    @root_validator
    def check_attributes_are_uniform(cls, values):
        if values["__root__"]:
            attrs = [tuple(s.data.keys()) for s in values["__root__"]]
            if len(set(attrs)) > 1:
                raise ValidationError("Attributes not all the same.")
            values["_attributes_str"] = ",".join(attrs[0])
        return values

    @classmethod
    def parse(cls, sol):
        pat = r"""
        :HRUStateVariableTable
        \s*:Attributes,(.+)
        \s*:Units.*?
        (.+)
        :EndHRUStateVariableTable
        """
        m = re.search(dedent(pat).strip(), sol, re.DOTALL)
        names = m.group(1).strip().split(",")
        lines = m.group(2).strip().splitlines()  # type: ignore

        states = [HRUState.parse(line, names) for line in lines]
        return cls.parse_obj(states)


class BasinIndex(Command):
    """Initial conditions for a flow segment."""

    # TODO: Check that qout and cie can be written with separating commas.
    index: int = 1
    name: str = "watershed"
    channel_storage: float = Field(0.0, alias="ChannelStorage")
    rivulet_storage: float = Field(0.0, alias="RivuletStorage")
    qout: Sequence[float] = Field((1.0, 0.0, 0.0), alias="Qout")
    qlat: Sequence[float] = Field(None, alias="Qlat")
    qin: Sequence[float] = Field(None, alias="Qin")

    _template = """
         :BasinIndex {index} {name}
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
        rec_values = {"index": index_name[0], "name": index_name[1]}
        for line in m.group(2).strip().splitlines():
            all_values = filter(None, re.split(r",|\s+", line.strip()))
            cmd, *values = all_values
            if cmd in [":ChannelStorage", ":RivuletStorage"]:
                assert len(values) == 1
                rec_values[cmd[1:]] = float(values[0])
            else:
                rec_values[cmd[1:]] = tuple(values)
        return cls(**rec_values)


class BasinStateVariables(Command):
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


class SoilClasses(Command):
    names: Sequence[str] = ()

    def commands(self):
        return "\n".join(self.names)


class VegetationClass(Command):
    name: str = ""
    max_ht: float = 0
    max_lai: float = 0
    max_leaf_cond: float = 0

    def to_rv(self):
        template = "{name:<16},{max_ht:>14},{max_lai:>14},{max_leaf_cond:>14}"
        return template.format(**self.dict())


class VegetationClasses(Command):
    __root__: Sequence[VegetationClass]
    _template = """
    :VegetationClasses
      :Attributes     ,        MAX_HT,       MAX_LAI, MAX_LEAF_COND
      :Units          ,             m,          none,      mm_per_s
    {_commands}
    :EndVegetationClasses
    """


class LandUseClass(Command):
    name: str = ""
    impermeable_frac: float = 0
    forest_coverage: float = 0

    def to_rv(self):
        template = "{name:<16},{impermeable_frac:>16},{forest_coverage:>16}"
        return template.format(**self.dict())

    def __str__(self):
        return self.to_rv()


class LandUseClasses(Command):
    __root__: Sequence[LandUseClass]
    _template = """
        :LandUseClasses
          :Attributes     ,IMPERMEABLE_FRAC, FOREST_COVERAGE
          :Units          ,           fract,           fract
        {_commands}
        :EndLandUseClasses
        """


@dataclass
class OldParameterList(Command):
    @dataclass
    class Record(Command):
        name: str = ""
        vals: Sequence[Union[float, None]] = ()

        @validator("vals", pre=True)
        def no_none_in_default(cls, v, values):
            """Make sure that no values are None for the [DEFAULT] record."""
            if values["name"] == "[DEFAULT]" and None in v:
                raise ValueError("Default record can not contain None.")
            return v

        def to_rv(self, **kwds):
            fmt = "{name:<16}" + len(self.vals) * ",{:>18}"
            evals = []
            for v in self.vals:
                ev = "_DEFAULT" if v is None else v
                evals.append(ev)

            return fmt.format(name=self.name, *evals)

    names: Sequence[str] = ()  # Subclass this with the right type.
    records: Sequence[Record] = ()

    _cmd: ClassVar[str]

    def to_rv(self, **kwds):
        template = """
        :{cmd}
            :Parameters     {parameter_names}
            :Units          {units}
            {records}
        :End{cmd}
        """

        fmt = ",{:>18}" * len(self.names)
        units = ",              none" * len(self.names)
        return dedent(template).format(
            cmd=self._cmd,
            parameter_names=fmt.format(*self.names),
            units=units,
            records="\n    ".join([r.to_rv(**kwds) for r in self.records]),
        )


class SoilParameterList(ParameterListCommand):
    names: Sequence[options.SoilParameters] = ()


class VegetationParameterList(ParameterListCommand):
    names: Sequence[options.VegetationParameters] = ()


class LandUseParameterList(ParameterListCommand):
    names: Sequence[options.LandUseParameters] = ()


# Aliases for convenience
# HRU = HRUsCommand.Record
# HRUState = HRUStateVariableTableCommand.Record
LU = LandUseClass
SB = SubBasin
SP = SoilProfile
VC = VegetationClasses
PL = ParameterList


# ---
