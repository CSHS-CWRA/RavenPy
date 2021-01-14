from dataclasses import dataclass, asdict, make_dataclass
from typing import List, NamedTuple, Tuple, Any


class Command:
    def __format__(self, format_spec):
        return self.to_rv()


@dataclass
class SubBasinsCommandRecord(Command):
    """Record to populate RVH :SubBasins command internal table."""

    subbasin_id: int = 0
    name: str = "sub_XXX"
    downstream_id: int = 0
    profile: str = "chn_XXX"
    reach_length: float = 0
    gauged: bool = False

    def to_rv(self):
        d = asdict(self)
        d["reach_length"] = d["reach_length"] if d["reach_length"] else "ZERO-"
        d["gauged"] = int(d["gauged"])
        return " ".join(map(str, d.values()))


@dataclass
class SubBasinsCommand(Command):
    """:SubBasins command (RVH)."""

    subbasins: Tuple[SubBasinsCommandRecord] = ()

    def to_rv(self):
        pat = """
:SubBasins
  :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
  :Units      none none          none    none           km    none
  {subbasin_records}
:EndSubBasins
        """
        recs = [f"\t{sb}" for sb in self.subbasins]
        return pat.format(subbasin_records="\n".join(recs))


@dataclass
class HRUsCommandRecord(Command):
    """Record to populate :HRUs command internal table (RVH)."""

    hru_id: int = 0
    area: float = 0  # km^2
    elevation: float = 0  # meters
    latitude: float = 0
    longitude: float = 0
    subbasin_id: int = 0
    land_use_class: str = ""
    veg_class: str = ""
    soil_profile: str = ""
    aquifer_profile: str = ""
    terrain_class: str = ""
    slope: float = 0.0
    aspect: float = 0.0

    def to_rv(self):
        d = asdict(self)
        return " ".join(map(str, d.values()))


@dataclass
class HRUsCommand(Command):
    """:HRUs command (RVH)."""

    hrus: Tuple[HRUsCommandRecord] = ()

    def to_rv(self):
        pat = """
:HRUs
  :Attributes      AREA  ELEVATION       LATITUDE      LONGITUDE BASIN_ID       LAND_USE_CLASS           VEG_CLASS      SOIL_PROFILE  AQUIFER_PROFILE TERRAIN_CLASS      SLOPE     ASPECT
  :Units            km2          m            deg            deg     none                  none               none              none             none          none        deg       degN
  {hru_records}
:EndHRUs
        """
        recs = [f"\t{hru}" for hru in self.hrus]
        return pat.format(hru_records="\n".join(recs))


@dataclass
class ReservoirCommand(Command):
    """:Reservoir command (RVH)."""

    subbasin_id: int = 0
    hru_id: int = 0
    name: str = "Lake_XXX"
    weir_coefficient: float = 0
    crest_width: float = 0
    max_depth: float = 0
    lake_area: float = 0

    def to_rv(self):
        pat = """
:Reservoir {name}
\t:SubBasinID {subbasin_id}
\t:HRUID {hru_id}
\t:Type RESROUTE_STANDARD
\t:WeirCoefficient {weir_coefficient}
\t:CrestWidth {crest_width}
\t:MaxDepth {max_depth}
\t:LakeArea {lake_area}
:EndReservoir
        """
        d = asdict(self)
        return pat.format(**d)


@dataclass
class ReservoirList(Command):
    """Sequence of :Reservoir commands."""

    reservoirs: Tuple[ReservoirCommand] = ()

    def to_rv(self):
        return "\n\n".join([r.to_rv() for r in self.reservoirs])


@dataclass
class SubBasinGroupCommand(Command):
    """:SubBasinGroup command (RVH)."""

    name: str = ""
    subbasin_ids: Tuple[int] = ()

    def to_rv(self):
        pat = """
:SubBasinGroup {name}
\t{subbasin_ids}
:EndSubBasinGroup
        """
        d = asdict(self)
        d["subbasin_ids"] = map(str, self.subbasin_ids)
        return pat.format(**d)


@dataclass
class ChannelProfileCommand(Command):
    """:ChannelProfile command (RVP)."""

    name: str = "chn_XXX"
    bed_slope: float = 0
    survey_points: Tuple[Tuple[float, float]] = ()
    roughness_zones: Tuple[Tuple[float, float]] = ()

    def to_rv(self):
        pat = """
:ChannelProfile {name}
\t:Bedslope {bed_slope}
\t:SurveyPoints
{survey_points}
\t:EndSurveyPoints
\t:RoughnessZones
{roughness_zones}
\t:EndRoughnessZones
:EndChannelProfile
        """
        d = asdict(self)
        d["survey_points"] = "\n".join(f"\t\t{p[0]} {p[1]}" for p in d["survey_points"])
        d["roughness_zones"] = "\n".join(
            f"\t\t{z[0]} {z[1]}" for z in d["roughness_zones"]
        )
        return pat.format(**d)


@dataclass
class ChannelProfileList(Command):
    channel_profiles: Tuple[ChannelProfileCommand] = ()

    def to_rv(self):
        return "\n\n".join([cp.to_rv() for cp in self.channel_profiles])


@dataclass
class GridWeightsCommand(Command):
    """:GridWeights command."""

    number_hrus: int = 0
    number_grid_cells: int = 0
    data: Tuple[Tuple[int, int, float]] = ()

    def to_rv(self, indent_level=0):
        indent = "\t" * indent_level
        pat = """
{indent}:GridWeights
{indent}\t:NumberHRUs {number_hrus}
{indent}\t:NumberGridCells {number_grid_cells}
{data}
{indent}:EndGridWeights
        """.strip()
        d = asdict(self)
        d["indent"] = indent
        d["data"] = "\n".join(f"{indent}\t{p[0]} {p[1]} {p[2]}" for p in self.data)
        return pat.format(**d)


@dataclass
class GriddedForcingCommand(Command):
    """:GriddedForcing command (RVT)."""

    name: str = ""
    forcing_type: str = ""
    file_name_nc: str = ""
    var_name_nc: str = ""
    dim_names_nc: Tuple[str, str, str] = ("x", "y", "t")
    grid_weights: GridWeightsCommand = None

    def to_rv(self):
        pat = """
:GriddedForcing {name}
\t:ForcingType {forcing_type}
\t:FileNameNC {file_name_nc}
\t:VarNameNC {var_name_nc}
\t:DimNamesNC {dim_names_nc_str}
{grid_weights}
:EndGriddedForcing
        """
        d = asdict(self)
        d["dim_names_nc_str"] = " ".join(self.dim_names_nc)
        d["grid_weights"] = self.grid_weights.to_rv()  # asdict seems to recurse... bummer
        return pat.format(**d)


@dataclass
class BaseValueCommand(Command):
    tag : str = ""
    value : Any = None

    def to_rv(self):
        pat = ":{tag} {value}"
        return pat.format(**asdict(self))


@dataclass
class BaseBooleanCommand(Command):
    tag: str = ""
    value: bool = False

    def to_rv(self):
        pat = ":{tag}"
        return pat.format(tag=self.tag) if self.value else ""


@dataclass
class LinearTransform(Command):
    scale: float = 1
    offset: float = 0

    def to_rv(self):
        pat = ":LinearTransform {slope:.15f} {intercept:.15f}"
        return pat.format(**asdict(self))


@dataclass
class DataCommand(Command):
    var = None
    path = None
    var_name = None
    units = None
    scale_factor = None
    add_offset = None
    time_shift = None
    dimensions = None
    index = None
    linear_transform = None

    runits = None
    raven_name = None
    site = None
    deaccumulate = None

    def to_rv(self):
        pat = """
:{kind} {raven_name} {site} {runits}
  :ReadFromNetCDF
     :FileNameNC      {path}
     :VarNameNC       {var_name}
     :DimNamesNC      {dimensions}
     :StationIdx      {index}
     {time_shift}
     {linear_transform}
     {deaccumulate}
  :EndReadFromNetCDF
:End{kind}"""
        d = asdict(self)
        d["kind"] = "ObservationData" if self.var == "water_volume_transport_in_river_channel" else "Data"


