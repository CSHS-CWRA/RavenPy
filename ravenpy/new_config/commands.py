import datetime as dt
import itertools
import re
from dataclasses import asdict, field
from itertools import chain
from pathlib import Path
from textwrap import dedent
from typing import (
    ClassVar,
    Dict,
    Literal,
    Optional,
    Sequence,
    Tuple,
    Union,
    no_type_check,
)

from pydantic import BaseModel, Field, HttpUrl, PrivateAttr, ValidationError, validator
from pydantic.dataclasses import dataclass

from ravenpy.config import options

from .base import Command, ParameterList, ParameterListCommand, RecordCommand

INDENT = " " * 4
VALUE_PADDING = 10


class SoilModel(Command):
    n: int = None

    def to_rv(self):
        return f":SoilModel            SOIL_MULTILAYER {self.n}\n"


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
    name: str
    soil_classes: Tuple[str, ...] = ()
    thicknesses: Tuple[float, ...] = ()

    def to_rv(self):
        # From the Raven manual: {profile_name,#horizons,{soil_class_name,thick.}x{#horizons}}x[NP]
        n_horizons = len(self.soil_classes)
        horizon_data = list(itertools.chain(*zip(self.soil_classes, self.thicknesses)))
        fmt = "{:<16},{:>4}," + ",".join(n_horizons * ["{:>12},{:>6}"])
        return fmt.format(self.name, n_horizons, *horizon_data)


class SoilProfiles(RecordCommand):
    record = SoilProfile


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


class SubBasins(RecordCommand):
    """SubBasins command (RVH)."""

    record = SubBasin
    template = """
        :SubBasins
          :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
          :Units      none none          none    none           km    none
          {records}
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


class HRUsCommand(RecordCommand):
    """HRUs command (RVH)."""

    record = HRU
    template = """
        :HRUs
          :Attributes      AREA  ELEVATION       LATITUDE      LONGITUDE BASIN_ID       LAND_USE_CLASS            VEG_CLASS      SOIL_PROFILE  AQUIFER_PROFILE TERRAIN_CLASS      SLOPE     ASPECT
          :Units            km2          m            deg            deg     none                  none                none              none             none          none        deg       degN
          {records}
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


@dataclass
class SBGroupPropertyMultiplierCommand(Command):

    group_name: str
    parameter_name: str
    mult: float

    def to_rv(self):
        template = ":SBGroupPropertyMultiplier {group_name} {parameter_name} {mult}"
        return dedent(template).format(**asdict(self))


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


class ReadFromNetCDF(Command):

    file_name_nc: Union[HttpUrl, Path] = Field(..., alias="FileNameNC")
    var_name_nc: str = Field(..., alias="VarNameNC")
    dim_names_nc: Sequence[str] = Field(..., alias="DimNamesNC")
    station_idx: int = Field(1, alias="StationIdx")
    time_shift: int = Field(None, alias="TimeShift")
    scale: float = 1
    offset: float = 0
    linear_transform: LinearTransform = Field(None, alias="LinearTransform")
    deaccumulate: bool = Field(None, alias="Deaccumulate")

    def __init__(self, **data):
        super().__init__(**data)
        if self.linear_transform is None:
            self.linear_transform = LinearTransform(
                scale=self.scale, offset=self.offset
            )

    @validator("dim_names_nc")
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


# class BaseData(Command):
#     """Do not use directly. Subclass."""
#
#     name: str = ""
#     units: str = ""
#     data_type: str = ""
#     # Can be a url or a path; note that the order "str, Path" is important here because
#     # otherwise pydantic would try to coerce a string into a Path, which is a problem for a url
#     # because it messes with slashes; so a Path will be one only when explicitly specified
#     file_name_nc: Union[HttpUrl, Path] = ""
#     var_name_nc: str = ""
#     dim_names_nc: Tuple[str, ...] = ("time",)
#     time_shift: RavenValue = Field(None, alias="TimeShift") # in days
#     scale: float = 1
#     offset: float = 0
#     deaccumulate: RavenBool = Field(None, alias="Deaccumulate")
#     latitude_var_name_nc: str = ""
#     longitude_var_name_nc: str = ""
#     elevation_var_name_nc: str = ""
#
#     @property
#     def dimensions(self):
#         """Return dimensions with time as the last dimension."""
#         dims = list(self.dim_names_nc)
#         for time_dim in ("t", "time"):
#             if time_dim in dims:
#                 dims.remove(time_dim)
#                 dims.append(time_dim)
#                 break
#         else:
#             raise Exception("No time dimension found in dim_names_nc tuple")
#         return " ".join(dims)
#
#     def asdict(self):
#         d = self.dict()
#         d["dimensions"] = self.dimensions
#         d["linear_transform"] = LinearTransform(scale=self.scale, offset=self.offset).to_rv()
#         d["deaccumulate"] = self.deaccumulate.to_rv("Deaccumulate")
#         d["time_shift"] = self.time_shift.to_rv("TimeShift")
#         if isinstance(d["file_name_nc"], Path):
#             # We can use the name of the file (as opposed to the full path)
#             # because we have a symlink to it in the execution folder
#             d["file_name_nc"] = d["file_name_nc"].name
#         return d


class Data(Command):
    data_type: str = ""
    site: str = ""
    units: str = ""
    read_from_netcdf: ReadFromNetCDF = Field(..., alias="ReadFromNetCDF")

    def to_rv(self):
        template = """
            :Data {data_type} {site} {units}
                :ReadFromNetCDF
                    :FileNameNC      {file_name_nc}
                    :VarNameNC       {var_name_nc}
                    :DimNamesNC      {dimensions}
                    :StationIdx      {index}
                    {time_shift}{linear_transform}{deaccumulate}
                :EndReadFromNetCDF
            :EndData
            """
        d = self.asdict()
        return dedent(template).format(**d)


class Gauge(Command):
    name: str = "default"
    latitude: float = 0
    longitude: float = 0
    elevation: float = 0

    # Accept strings to embed parameter names into Ostrich templates
    rain_correction: float = Field(None, alias="RainCorrection")
    snow_correction: float = Field(None, alias="SnowCorrection")

    monthly_ave_evaporation: Tuple[float, ...] = ()
    monthly_ave_temperature: Tuple[float, ...] = ()

    data_cmds: Tuple[Data, ...] = ()

    def to_rv(self):
        template = """
        :Gauge {name}
          :Latitude {latitude}
          :Longitude {longitude}
          :Elevation {elevation}
          {rain_correction}{snow_correction}{monthly_ave_evaporation}{monthly_ave_temperature}
          {data_cmds}
        :EndGauge
        """

        d = self.dict()
        d["rain_correction"] = (
            f":RainCorrection {self.rain_correction}\n" if self.rain_correction else ""
        )
        d["snow_correction"] = (
            f":SnowCorrection {self.snow_correction}\n" if self.snow_correction else ""
        )
        if self.monthly_ave_evaporation:
            evap_data = " ".join(map(str, self.monthly_ave_evaporation))
            d["monthly_ave_evaporation"] = f":MonthlyAveEvaporation {evap_data}\n"
        else:
            d["monthly_ave_evaporation"] = ""
        if self.monthly_ave_temperature:
            temp_data = " ".join(map(str, self.monthly_ave_temperature))
            d["monthly_ave_temperature"] = f":MonthlyAveTemperature {temp_data}\n"
        else:
            d["monthly_ave_temperature"] = ""
        d["data_cmds"] = "\n\n".join(map(str, self.data_cmds))  # type: ignore
        return dedent(template).format(**d)


class ObservationDataCommand(Data):
    subbasin_id: int = 1

    def to_rv(self):
        template = """
        :ObservationData {data_type} {subbasin_id} {units}
            :ReadFromNetCDF
                :FileNameNC      {file_name_nc}
                :VarNameNC       {var_name_nc}
                :DimNamesNC      {dimensions}
                :StationIdx      {index}
                {time_shift}{linear_transform}{deaccumulate}
            :EndReadFromNetCDF
        :EndObservationData
        """
        d = self.asdict()
        return dedent(template).format(**d)


class GridWeightsCommand(Command):
    """GridWeights command.

    Important note: this command can be embedded in both a `GriddedForcingCommand`
    or a `StationForcingCommand`.

    The default is to have a single cell that covers an entire single HRU, with a
    weight of 1.
    """

    number_hrus: int = 1
    number_grid_cells: int = 1
    data: Tuple[Tuple[int, int, float], ...] = ((1, 0, 1.0),)

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

    def to_rv(self, indent_level=0):
        template = """
        {indent}:GridWeights
        {indent}    :NumberHRUs {number_hrus}
        {indent}    :NumberGridCells {number_grid_cells}
        {data}
        {indent}:EndGridWeights
        """
        indent = INDENT * indent_level
        d = asdict(self)
        d["indent"] = indent
        d["data"] = "\n".join(f"{indent}    {p[0]} {p[1]} {p[2]}" for p in self.data)
        return dedent(template).strip().format(**d)


class RedirectToFileCommand(Command):
    """RedirectToFile command (RVT).

    For the moment, this command can only be used in the context of a `GriddedForcingCommand`
    or a `StationForcingCommand`, as a `grid_weights` field replacement when inlining is not
    desired.
    """

    path: Path

    def to_rv(self, indent_level=0):
        template = "{indent}:RedirectToFile {path}"
        indent = INDENT * indent_level
        d = asdict(self)
        d["indent"] = indent
        # We can use the name of the file (as opposed to the full path)
        # because we have a symlink to it in the execution folder
        d["path"] = d["path"].name
        return template.format(**d)


class GriddedForcingCommand(BaseData):
    """GriddedForcing command (RVT)."""

    dim_names_nc: Tuple[str, str, str] = ("x", "y", "t")
    grid_weights: Union[
        GridWeightsCommand, RedirectToFileCommand
    ] = GridWeightsCommand()

    # :LatitudeVarNameNC {latitude_var_name_nc}
    # :LongitudeVarNameNC {longitude_var_name_nc}
    # :ElevationVarNameNC {elevation_var_name_nc}

    def to_rv(self):
        template = """
            :GriddedForcing {name}
                :ForcingType {data_type}
                :FileNameNC {file_name_nc}
                :VarNameNC {var_name_nc}
                :DimNamesNC {dimensions}
                {time_shift}{linear_transform}{deaccumulate}
            {grid_weights}
            :EndGriddedForcing
            """
        d = self.asdict()
        d["grid_weights"] = self.grid_weights.to_rv(indent_level=1)
        return dedent(template).format(**d)


class StationForcingCommand(BaseData):
    """StationForcing command (RVT)."""

    dim_names_nc: Tuple[str, str] = ("station", "time")
    grid_weights: Union[
        GridWeightsCommand, RedirectToFileCommand
    ] = GridWeightsCommand()

    # :LatitudeVarNameNC {latitude_var_name_nc}
    # :LongitudeVarNameNC {longitude_var_name_nc}
    # :ElevationVarNameNC {elevation_var_name_nc}

    def to_rv(self):
        template = """
            :StationForcing {name} {units}
                :ForcingType {data_type}
                :FileNameNC {file_name_nc}
                :VarNameNC {var_name_nc}
                :DimNamesNC {dimensions}
                {time_shift}{linear_transform}{deaccumulate}
            {grid_weights}
            :EndStationForcing
            """
        d = self.asdict()
        d["grid_weights"] = self.grid_weights.to_rv(indent_level=1)
        return dedent(template).format(**d)


# class HRUStateVariableTableCommand(Command):
#     """Initial condition for a given HRU."""
#
#     class Record(Command):
#         index: int = 1
#         data: Dict[str, Union[float, str]] = field(default_factory=dict)
#
#         def to_rv(self):
#             return ",".join(map(str, (self.index,) + tuple(self.data.values())))
#
#     hru_states: Dict[int, Record] = field(default_factory=dict)
#
#     @classmethod
#     def parse(cls, sol):
#         pat = r"""
#         :HRUStateVariableTable
#         \s*:Attributes,(.+)
#         \s*:Units.*?
#         (.+)
#         :EndHRUStateVariableTable
#         """
#         m = re.search(dedent(pat).strip(), sol, re.DOTALL)
#         names = m.group(1).strip().split(",")
#         lines = m.group(2).strip().splitlines()  # type: ignore
#         lines = [re.split(r",|\s+", line.strip()) for line in lines]
#         hru_states = {}
#         for line in lines:
#             idx, *values = line
#             idx = int(idx)
#             values = list(map(float, values))
#             hru_states[idx] = cls.Record(index=idx, data=dict(zip(names, values)))
#         return cls(hru_states)
#
#     def to_rv(self):
#         template = """
#             :HRUStateVariableTable
#                 :Attributes,{names}
#                 {values}
#             :EndHRUStateVariableTable
#             """
#         names = sorted(
#             list(set(chain(*[tuple(s.data.keys()) for s in self.hru_states.values()])))
#         )
#         values = [
#             [
#                 s.index,
#             ]
#             + [s.data.get(n, 0.0) for n in names]
#             for s in self.hru_states.values()
#         ]
#         return dedent(template).format(
#             names=",".join(names),
#             values="\n    ".join([",".join(map(str, v)) for v in values]),
#         )
#
#
# class BasinIndexCommand(Command):
#     """Initial conditions for a flow segment."""
#
#     index: int = 1
#     name: str = "watershed"
#     channel_storage: float = 0
#     rivulet_storage: float = 0
#     qout: Tuple[float, ...] = (1, 0, 0)
#     qin: Optional[Tuple[float, ...]] = None
#     qlat: Optional[Tuple[float, ...]] = None
#
#     @classmethod
#     @no_type_check
#     def parse(cls, s):
#         pat = r"""
#         :BasinIndex (.+?)
#         (.+)
#         """
#         m = re.search(dedent(pat).strip(), s, re.DOTALL)
#         index_name = re.split(r",|\s+", m.group(1).strip())
#         rec_values = {"index": index_name[0], "name": index_name[1]}
#         for line in m.group(2).strip().splitlines():
#             all_values = filter(None, re.split(r",|\s+", line.strip()))
#             cmd, *values = all_values
#             if cmd == ":ChannelStorage":
#                 assert len(values) == 1
#                 rec_values["channel_storage"] = float(values[0])
#             elif cmd == ":RivuletStorage":
#                 assert len(values) == 1
#                 rec_values["rivulet_storage"] = float(values[0])
#             else:
#                 rec_values[cmd[1:].lower()] = tuple(values)
#         return cls(**rec_values)
#
#     def to_rv(self):
#         template = """
#         :BasinIndex {index} {name}
#             :ChannelStorage {channel_storage}
#             :RivuletStorage {rivulet_storage}
#             {qout}
#             {qin}
#             {qlat}
#             """
#         d = asdict(self)
#         for k in ["qout", "qin", "qlat"]:
#             if d[k]:
#                 v = " ".join(map(str, d[k]))
#                 q = k.capitalize()
#                 d[k] = f":{q} {v}"
#             else:
#                 d[k] = ""
#         return dedent(template).format(**d)
#
#
# class BasinStateVariablesCommand(Command):
#
#     basin_states: Dict[int, BasinIndexCommand] = field(default_factory=dict)
#
#     @classmethod
#     @no_type_check
#     def parse(cls, sol):
#         pat = r"""
#         :BasinStateVariables
#         (.+)
#         :EndBasinStateVariables
#         """
#         m = re.search(dedent(pat).strip(), sol, re.DOTALL)
#         bi_strings = filter(None, m.group(1).strip().split(":BasinIndex"))
#         basin_states = {}
#         for bi_string in bi_strings:
#             bi = BasinIndexCommand.parse(f":BasinIndex {bi_string}")
#             basin_states[bi.index] = bi
#         return cls(basin_states)
#
#     def to_rv(self):
#         template = """
#             :BasinStateVariables
#                 {basin_states_list}
#             :EndBasinStateVariables
#             """
#
#         return dedent(template).format(
#             basin_states_list="\n".join(map(str, self.basin_states.values()))
#         )
#


class SoilClasses(RecordCommand):
    record: str


class VegetationClasses(Command):
    name: str = ""
    max_ht: float = 0
    max_lai: float = 0
    max_leaf_cond: float = 0

    def to_rv(self):
        template = "{name:<16},{max_ht:>14},{max_lai:>14},{max_leaf_cond:>14}"
        return template.format(**asdict(self))


class VegetationClassesCommand(RecordCommand):
    record = VegetationClasses
    template = """
    :VegetationClasses
        :Attributes     ,        MAX_HT,       MAX_LAI, MAX_LEAF_COND
        :Units          ,             m,          none,      mm_per_s
        {records}
    :EndVegetationClasses
    """


class LandUseClass(Command):
    name: str = ""
    impermeable_frac: float = 0
    forest_coverage: float = 0

    def to_rv(self):
        template = "{name:<16},{impermeable_frac:>16},{forest_coverage:>16}"
        return template.format(**asdict(self))

    def __str__(self):
        return self.to_rv()


class LandUseClassesCommand(RecordCommand):
    record = LandUseClass
    template = """
        :LandUseClasses
            :Attributes     ,IMPERMEABLE_FRAC, FOREST_COVERAGE
            :Units          ,           fract,           fract
            {records}
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
    records: Sequence[ParameterList] = ()


# Aliases for convenience
# HRU = HRUsCommand.Record
# HRUState = HRUStateVariableTableCommand.Record
LU = LandUseClass
SB = SubBasin
SP = SoilProfile
VC = VegetationClasses
PL = ParameterList


# ---
