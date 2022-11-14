import datetime as dt
import itertools
import re
from dataclasses import asdict, field
from itertools import chain
from pathlib import Path
from textwrap import dedent
from typing import ClassVar, Dict, Optional, Sequence, Tuple, Union, no_type_check

from pydantic import validator
from pydantic.dataclasses import dataclass

from . import options
from .base import (
    RavenCoefficient,
    RavenCommand,
    RavenOption,
    RavenOptionList,
    RavenSwitch,
    RavenValue,
)

INDENT = " " * 4
VALUE_PADDING = 10


# --- Boolean switches --- #
class DebugMode(RavenSwitch):
    """"""


class DeltaresFEWSMode(RavenSwitch):
    """"""


class DirectEvaporation(RavenSwitch):
    """Rainfall is automatically reduced through evapotranspiration up to the limit of the calculated PET."""


class DontWriteWatershedStorage(RavenSwitch):
    """Do not write watershed storage variables to disk."""


class NetCDFAttribute(RavenSwitch):
    """"""


class NoisyMode(RavenSwitch):
    """"""


class PavicsMode(RavenSwitch):
    """"""


class SilentMode(RavenSwitch):
    """"""


class SuppressOutput(RavenSwitch):
    """Write minimal output to disk when enabled."""


class WriteForcingFunctions(RavenSwitch):
    """Write watershed averaged forcing functions (e.g. rainfall, radiation, PET, etc)."""


class WriteSubbasinFile(RavenSwitch):
    """"""


# --- Coefficients --- #


class AirSnowCoeff(RavenCoefficient):
    """The air/snow heat transfer coefficient as used in the `SNOTEMP_NEWTONS` snow temperature evolution routine."""

    value: float
    """Heat transfer coefficient [1/d]."""


class AvgAnnualSnow(RavenCoefficient):
    """The average annual snow for the entire watershed used in the CEMANEIGE algorithm."""

    value: float
    """Average annual snow [mm]."""


class PrecipitationLapseRate(RavenCoefficient):
    """The simple linear precipitation lapse rate  used in the `OROCORR_SIMPLELAPSE` orographic correction algorithm."""

    value: float
    """Lapse rate [mm/d/km]"""


class AdiabaticLapseRate(RavenCoefficient):
    """Base adiabatic lapse rate."""

    value: float
    """Base adiabatic lapse rate [C/km]"""


# --- Options --- #


@dataclass
class Calendar(RavenOption):
    option: options.Calendar = options.Calendar.STANDARD


@dataclass
class CatchmentRoute(RavenOption):
    option: options.CatchmentRoute


@dataclass
class CloudCoverMethod(RavenOption):
    option: options.CloudCoverMethod


@dataclass
class Duration(RavenValue):
    """Duration of simulation.

    Attributes
    ----------
    values : float
        Simulation duration [d].
    """


@dataclass
class EndDate(RavenValue):
    value: dt.datetime = None


@dataclass
class EvaluationMetrics(RavenOptionList):
    options: Sequence[options.EvaluationMetrics]


@dataclass
class Evaporation(RavenOption):
    option: options.Evaporation


@dataclass
class LakeStorage(RavenValue):
    value: options.StateVariables


@dataclass
class OW_Evaporation(RavenOption):
    option: options.Evaporation


@dataclass
class MonthlyInterpolationMethod(RavenOption):
    option: options.MonthlyInterpolationMethod


@dataclass
class PotentialMeltMethod(RavenOption):
    """Potential snow melt algorithm."""

    option: options.PotentialMeltMethod
    """Potential melt algorithm."""


@dataclass
class PrecipIceptFract(RavenOption):
    option: options.PrecipIceptFract


@dataclass
class RainSnowFraction(RavenOption):
    option: options.RainSnowFraction = options.RainSnowFraction.DATA


@dataclass
class RelativeHumidityMethod(RavenOption):
    option: options.RelativeHumidityMethod


@dataclass
class Routing(RavenOption):
    option: options.Routing = options.Routing.NONE


@dataclass
class RunName(RavenValue):
    value: str = "run"


@dataclass
class SoilModel(RavenOption):
    option: options.SoilModel
    n: int = None

    def to_rv(self):
        n = self.n if self.n is not None else ""
        return f":SoilModel           {self.options.value} {n}\n"


@dataclass
class StartDate(RavenValue):
    value: dt.datetime = None


@dataclass
class SubdailyMethod(RavenOption):
    option: options.SubdailyMethod


@dataclass
class SWCanopyCorrect(RavenOption):
    option: options.SWCanopyCorrect


@dataclass
class SWCloudCorrect(RavenOption):
    option: options.SWCloudCorrect


@dataclass
class SWRadiationMethod(RavenOption):
    option: options.SWRadiationMethod


@dataclass
class TimeStep(RavenValue):
    value: float = 1.0


@dataclass
class WindspeedMethod(RavenOption):
    option: options.WindspeedMethod


# --- Custom commands --- #


@dataclass
class LinearTransform(RavenCommand):
    scale: Optional[float] = 1
    offset: Optional[float] = 0

    def to_rv(self):
        template = ":LinearTransform {scale:.15f} {offset:.15f}\n"
        if (self.scale != 1) or (self.offset != 0):
            return template.format(**asdict(self))
        return ""


@dataclass
class RainSnowTransition(RavenCommand):
    """Specify the range of temperatures over which there will be a rain/snow mix when partitioning total precipitation into rain and snow components."""

    temp: float
    """Midpoint of the temperature range [C]."""
    delta: float
    """Range [C]."""

    def to_rv(self):
        template = ":RainSnowTransition {temp} {float}"
        return template.format(**asdict(self))


@dataclass
class EvaluationPeriod(RavenCommand):
    """:EvaluationPeriod [period_name] [start yyyy-mm-dd] [end yyyy-mm-dd]"""

    name: str
    start: dt.date
    end: dt.date

    def to_rv(self):
        template = ":EvaluationPeriod {name} {start} {end}"
        return template.format(**asdict(self))


@dataclass
class CustomOutput(RavenCommand):
    """Create custom output file to track a single variable, parameter or forcing function over time at a number of basins, HRUs, or across the watershed."""

    time_per: str
    """Time period. Allowed keys: {'DAILY', 'MONTHLY', 'YEARLY', 'WATER_YEARLY', 'CONTINUOUS'}"""
    stat: str
    """
    Statistic reported for each time interval.
    Allowed keys: {'AVERAGE', 'MAXIMUM', 'MINIMUM', 'RANGE', 'MEDIAN', 'QUARTILES', 'HISTOGRAM [min] [max] [# bins]'}
    """
    variable: str
    """Variable or parameter name. Consult the Raven documentation for the list of allowed names."""
    space_agg: str
    """
    Spatial evaluation domain.
    Allowed keys: {'BY_BASIN', 'BY_HRU', 'BY_HRU_GROUP', 'BY_SB_GROUP', 'ENTIRE_WATERSHED'}
    """
    filename: str = ""
    """
    Output file name.
    Defaults to something approximately like `<run name>_<variable>_<time_per>_<stat>_<space_agg>.nc`
    """

    def to_rv(self):
        template = ":CustomOutput {time_per} {stat} {variable} {space_agg} {filename}"
        return template.format(**asdict(self))


@dataclass
class SubBasinsCommand(RavenCommand):
    """SubBasins command (RVH)."""

    @dataclass
    class Record(RavenCommand):
        """Record to populate RVH :SubBasins command internal table."""

        subbasin_id: int = 0
        name: str = "sub_XXX"
        downstream_id: int = 0
        profile: str = "chn_XXX"
        reach_length: float = 0
        gauged: bool = False
        gauge_id: Optional[str] = ""  # This attribute is not rendered to RVH

        def to_rv(self):
            d = asdict(self)
            d["reach_length"] = d["reach_length"] if d["reach_length"] else "ZERO-"
            d["gauged"] = int(d["gauged"])
            del d["gauge_id"]
            return " ".join(f"{v: <{VALUE_PADDING}}" for v in d.values())

    subbasins: Tuple[Record, ...] = ()

    def to_rv(self):
        template = """
            :SubBasins
                :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
                :Units      none none          none    none           km    none
            {subbasin_records}
            :EndSubBasins
        """
        recs = [f"    {sb}" for sb in self.subbasins]
        return dedent(template).format(subbasin_records="\n".join(recs))


@dataclass
class HRUsCommand(RavenCommand):
    """HRUs command (RVH)."""

    @dataclass
    class Record(RavenCommand):
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
            d = asdict(self)
            del d["hru_type"]
            return " ".join(f"{v: <{VALUE_PADDING * 2}}" for v in d.values())

    hrus: Tuple[Record, ...] = ()

    def to_rv(self):
        template = """
            :HRUs
                :Attributes      AREA  ELEVATION       LATITUDE      LONGITUDE BASIN_ID       LAND_USE_CLASS            VEG_CLASS      SOIL_PROFILE  AQUIFER_PROFILE TERRAIN_CLASS      SLOPE     ASPECT
                :Units            km2          m            deg            deg     none                  none                none              none             none          none        deg       degN
            {hru_records}
            :EndHRUs
            """
        recs = [f"    {hru}" for hru in self.hrus]
        return dedent(template).format(hru_records="\n".join(recs))


@dataclass
class ReservoirCommand(RavenCommand):
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
        d = asdict(self)
        return dedent(template).format(**d)


@dataclass
class SubBasinGroupCommand(RavenCommand):
    """SubBasinGroup command (RVH)."""

    name: str = ""
    subbasin_ids: Tuple[int, ...] = ()

    def to_rv(self):
        template = """
            :SubBasinGroup {name}
                {subbasin_ids}
            :EndSubBasinGroup
            """
        d = asdict(self)
        n_per_line = 10
        sbids = sorted(self.subbasin_ids)
        sbids_lines = [
            map(str, sbids[i : i + n_per_line])
            for i in range(0, len(sbids), n_per_line)
        ]
        d["subbasin_ids"] = "\n    ".join([" ".join(sbids) for sbids in sbids_lines])
        return dedent(template).format(**d)


@dataclass
class SBGroupPropertyMultiplierCommand(RavenCommand):
    """"""

    group_name: str
    parameter_name: str
    mult: float

    def to_rv(self):
        template = ":SBGroupPropertyMultiplier {group_name} {parameter_name} {mult}"
        return dedent(template).format(**asdict(self))


@dataclass
class ChannelProfileCommand(RavenCommand):
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
        d = asdict(self)
        d["survey_points"] = "\n".join(
            f"{INDENT * 2}{p[0]} {p[1]}" for p in d["survey_points"]
        )
        d["roughness_zones"] = "\n".join(
            f"{INDENT * 2}{z[0]} {z[1]}" for z in d["roughness_zones"]
        )
        return dedent(template).format(**d)


@dataclass
class BaseDataCommand(RavenCommand):
    """Do not use directly. Subclass."""

    name: Optional[str] = ""
    units: Optional[str] = ""
    data_type: str = ""
    # Can be an url or a path; note that the order "str, Path" is important here because
    # otherwise pydantic would try to coerce a string into a Path, which is a problem for a url
    # because it messes with slashes; so a Path will be one only when explicitly specified
    file_name_nc: Union[str, Path] = ""  # can be a URL
    var_name_nc: str = ""
    dim_names_nc: Tuple[str, ...] = ("time",)
    time_shift: Optional[float] = None  # in days
    scale: float = 1
    offset: float = 0
    deaccumulate: Optional[bool] = False
    latitude_var_name_nc: str = ""
    longitude_var_name_nc: str = ""
    elevation_var_name_nc: str = ""

    @property
    def dimensions(self):
        """Return dimensions with time as the last dimension."""
        dims = list(self.dim_names_nc)
        for time_dim in ("t", "time"):
            if time_dim in dims:
                dims.remove(time_dim)
                dims.append(time_dim)
                break
        else:
            raise Exception("No time dimension found in dim_names_nc tuple")
        return " ".join(dims)

    @property
    def linear_transform(self):
        return LinearTransform(scale=self.scale, offset=self.offset)

    def asdict(self):
        d = asdict(self)
        d["dimensions"] = self.dimensions
        d["linear_transform"] = self.linear_transform
        d["deaccumulate"] = ":Deaccumulate\n" if self.deaccumulate else ""
        d["time_shift"] = f":TimeShift {self.time_shift}\n" if self.time_shift else ""
        if isinstance(d["file_name_nc"], Path):
            # We can use the name of the file (as opposed to the full path)
            # because we have a symlink to it in the execution folder
            d["file_name_nc"] = d["file_name_nc"].name
        return d


@dataclass
class DataCommand(BaseDataCommand):
    index: int = 1  # Indexing starts with 1.
    data_type: str = ""
    site: str = ""
    var: str = ""

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


@dataclass
class GaugeCommand(RavenCommand):
    name: str = "default"
    latitude: float = 0
    longitude: float = 0
    elevation: float = 0

    # Accept strings to embed parameter names into Ostrich templates
    rain_correction: Optional[Union[float, str]] = 1
    snow_correction: Optional[Union[float, str]] = 1

    monthly_ave_evaporation: Optional[Tuple[float, ...]] = ()
    monthly_ave_temperature: Optional[Tuple[float, ...]] = ()

    data_cmds: Optional[Tuple[DataCommand, ...]] = ()

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

        d = asdict(self)
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


@dataclass
class ObservationDataCommand(DataCommand):
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


@dataclass
class GridWeightsCommand(RavenCommand):
    """GridWeights command.

    Important note: this command can be embedded in both a `GriddedForcingCommand` or a `StationForcingCommand`.
    The default is to have a single cell that covers an entire single HRU, with a weight of 1.
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


@dataclass
class RedirectToFileCommand(RavenCommand):
    """RedirectToFile command (RVT).

    For the moment, this command can only be used in the context of a `GriddedForcingCommand` or a
    `StationForcingCommand`, as a `grid_weights` field replacement when inlining is not desired.
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


@dataclass
class GriddedForcingCommand(BaseDataCommand):
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


@dataclass
class StationForcingCommand(BaseDataCommand):
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


@dataclass
class HRUStateVariableTableCommand(RavenCommand):
    """Initial condition for a given HRU."""

    @dataclass
    class Record(RavenCommand):
        index: int = 1
        data: Dict[str, Union[float, str]] = field(default_factory=dict)

        def to_rv(self):
            return ",".join(map(str, (self.index,) + tuple(self.data.values())))

    hru_states: Dict[int, Record] = field(default_factory=dict)

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
        lines = [re.split(r",|\s+", line.strip()) for line in lines]
        hru_states = {}
        for line in lines:
            idx, *values = line
            idx = int(idx)
            values = list(map(float, values))
            hru_states[idx] = cls.Record(index=idx, data=dict(zip(names, values)))
        return cls(hru_states)

    def to_rv(self):
        template = """
            :HRUStateVariableTable
                :Attributes,{names}
                {values}
            :EndHRUStateVariableTable
            """
        names = sorted(
            list(set(chain(*[tuple(s.data.keys()) for s in self.hru_states.values()])))
        )
        values = [
            [
                s.index,
            ]
            + [s.data.get(n, 0.0) for n in names]
            for s in self.hru_states.values()
        ]
        return dedent(template).format(
            names=",".join(names),
            values="\n    ".join([",".join(map(str, v)) for v in values]),
        )


@dataclass
class BasinIndexCommand(RavenCommand):
    """Initial conditions for a flow segment."""

    index: int = 1
    name: str = "watershed"
    channel_storage: float = 0
    rivulet_storage: float = 0
    qout: Tuple[float, ...] = (1, 0, 0)
    qin: Optional[Tuple[float, ...]] = None
    qlat: Optional[Tuple[float, ...]] = None

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
            if cmd == ":ChannelStorage":
                assert len(values) == 1
                rec_values["channel_storage"] = float(values[0])
            elif cmd == ":RivuletStorage":
                assert len(values) == 1
                rec_values["rivulet_storage"] = float(values[0])
            else:
                rec_values[cmd[1:].lower()] = tuple(values)
        return cls(**rec_values)

    def to_rv(self):
        template = """
        :BasinIndex {index} {name}
            :ChannelStorage {channel_storage}
            :RivuletStorage {rivulet_storage}
            {qout}
            {qin}
            {qlat}
            """
        d = asdict(self)
        for k in ["qout", "qin", "qlat"]:
            if d[k]:
                v = " ".join(map(str, d[k]))
                q = k.capitalize()
                d[k] = f":{q} {v}"
            else:
                d[k] = ""
        return dedent(template).format(**d)


@dataclass
class BasinStateVariablesCommand(RavenCommand):

    basin_states: Dict[int, BasinIndexCommand] = field(default_factory=dict)

    @classmethod
    @no_type_check
    def parse(cls, sol):
        pat = r"""
        :BasinStateVariables
        (.+)
        :EndBasinStateVariables
        """
        m = re.search(dedent(pat).strip(), sol, re.DOTALL)
        bi_strings = filter(None, m.group(1).strip().split(":BasinIndex"))
        basin_states = {}
        for bi_string in bi_strings:
            bi = BasinIndexCommand.parse(f":BasinIndex {bi_string}")
            basin_states[bi.index] = bi
        return cls(basin_states)

    def to_rv(self):
        template = """
            :BasinStateVariables
                {basin_states_list}
            :EndBasinStateVariables
            """

        return dedent(template).format(
            basin_states_list="\n".join(map(str, self.basin_states.values()))
        )


@dataclass
class SoilClassesCommand(RavenCommand):
    @dataclass
    class Record(RavenCommand):
        name: str = ""

        def to_rv(self):
            return " ".join(map(str, asdict(self).values()))

    soil_classes: Tuple[Record, ...] = ()

    def to_rv(self):
        template = """
            :SoilClasses
                {soil_class_records}
            :EndSoilClasses
            """
        return dedent(template).format(
            soil_class_records="\n    ".join(map(str, self.soil_classes))
        )


@dataclass
class SoilProfilesCommand(RavenCommand):
    @dataclass
    class Record(RavenCommand):
        profile_name: str = ""
        soil_class_names: Tuple[str, ...] = ()
        thicknesses: Tuple[float, ...] = ()

        def to_rv(self):
            # From the Raven manual: {profile_name,#horizons,{soil_class_name,thick.}x{#horizons}}x[NP]
            n_horizons = len(self.soil_class_names)
            horizon_data = list(
                itertools.chain(*zip(self.soil_class_names, self.thicknesses))
            )
            fmt = "{:<16},{:>4}," + ",".join(n_horizons * ["{:>12},{:>6}"])
            return fmt.format(self.profile_name, n_horizons, *horizon_data)

    soil_profiles: Tuple[Record, ...] = ()

    def to_rv(self):
        template = """
            :SoilProfiles
                {soil_profile_records}
            :EndSoilProfiles
            """
        return dedent(template).format(
            soil_profile_records="\n".join(map(str, self.soil_profiles))
        )


@dataclass
class VegetationClassesCommand(RavenCommand):
    @dataclass
    class Record(RavenCommand):
        name: str = ""
        max_ht: float = 0
        max_lai: float = 0
        max_leaf_cond: float = 0

        def to_rv(self):
            template = "{name:<16},{max_ht:>14},{max_lai:>14},{max_leaf_cond:>14}"
            return template.format(**asdict(self))

    vegetation_classes: Tuple[Record, ...] = ()

    def to_rv(self):
        template = """
        :VegetationClasses
            :Attributes     ,        MAX_HT,       MAX_LAI, MAX_LEAF_COND
            :Units          ,             m,          none,      mm_per_s
            {vegetation_class_records}
        :EndVegetationClasses
        """
        return dedent(template).format(
            vegetation_class_records="\n".join(map(str, self.vegetation_classes))
        )


@dataclass
class LandUseClassesCommand(RavenCommand):
    @dataclass
    class Record(RavenCommand):
        name: str = ""
        impermeable_frac: float = 0
        forest_coverage: float = 0

        def to_rv(self):
            template = "{name:<16},{impermeable_frac:>16},{forest_coverage:>16}"
            return template.format(**asdict(self))

    land_use_classes: Tuple[Record, ...] = ()

    def to_rv(self):
        template = """
            :LandUseClasses
                :Attributes     ,IMPERMEABLE_FRAC, FOREST_COVERAGE
                :Units          ,           fract,           fract
                {land_use_class_records}
            :EndLandUseClasses
            """
        return dedent(template).format(
            land_use_class_records="\n".join(map(str, self.land_use_classes))
        )


@dataclass
class ParameterList(RavenCommand):
    @dataclass
    class Record(RavenCommand):
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


@dataclass
class SoilParameterListCommand(ParameterList):
    names: Sequence[options.SoilParameters] = ()
    _cmd = "SoilParameterList"


@dataclass
class VegetationParameterListCommand(ParameterList):
    names: Sequence[options.VegetationParameters] = ()
    _cmd = "VegetationParameterList"


@dataclass
class LandUseParameterListCommand(ParameterList):
    names: Sequence[options.LandUseParameters] = ()
    _cmd = "LandUseParameterList"


# Aliases for convenience
HRU = HRUsCommand.Record
HRUState = HRUStateVariableTableCommand.Record
LU = LandUseClassesCommand.Record
Sub = SubBasinsCommand.Record
SOIL = SoilProfilesCommand.Record
VEG = VegetationClassesCommand.Record
PL = ParameterList.Record
