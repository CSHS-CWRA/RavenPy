from dataclasses import dataclass
from typing import List, NamedTuple, Tuple


class SubBasinsCommandRecord(NamedTuple):
    """Record to populate RVH :SubBasins command internal table."""

    subbasin_id: int = 0
    name: str = "sub_XXX"
    downstream_id: int = 0
    profile: str = "chn_XXX"
    reach_length: float = 0
    gauged: bool = False

    def to_rv(self):
        d = self._asdict()
        d["reach_length"] = d["reach_length"] if d["reach_length"] else "ZERO-"
        d["gauged"] = int(d["gauged"])
        return " ".join(map(str, d.values()))


class SubBasinsCommand(NamedTuple):
    """:SubBasins command (RVH)."""

    subbasins: List[SubBasinsCommandRecord] = []

    def to_rv(self):
        pat = """
:SubBasins
  :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
  :Units      none none          none    none           km    none
  {subbasin_records}
:EndSubBasins
        """
        recs = [f"\t{sb.to_rv()}" for sb in self.subbasins]
        return pat.format(subbasin_records="\n".join(recs))


class HRUsCommandRecord(NamedTuple):
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
        d = self._asdict()
        return " ".join(map(str, d.values()))


class HRUsCommand(NamedTuple):
    """:HRUs command (RVH)."""

    hrus: List[HRUsCommandRecord] = []

    def to_rv(self):
        pat = """
:HRUs
  :Attributes      AREA  ELEVATION       LATITUDE      LONGITUDE BASIN_ID       LAND_USE_CLASS           VEG_CLASS      SOIL_PROFILE  AQUIFER_PROFILE TERRAIN_CLASS      SLOPE     ASPECT
  :Units            km2          m            deg            deg     none                  none               none              none             none          none        deg       degN
  {hru_records}
:EndHRUs
        """
        recs = [f"\t{hru.to_rv()}" for hru in self.hrus]
        return pat.format(hru_records="\n".join(recs))


class ReservoirCommand(NamedTuple):
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
        d = self._asdict()
        return pat.format(**d)


class SubBasinGroupCommand(NamedTuple):
    """:SubBasinGroup command (RVH)."""

    name: str = ""
    subbasin_ids: List[int] = []

    def to_rv(self):
        pat = """
:SubBasinGroup {name}
\t{subbasin_ids}
:EndSubBasinGroup
        """
        d = self._asdict()
        d["subbasin_ids"] = map(str, self.subbasin_ids)
        return pat.format(**d)


class ChannelProfileCommand(NamedTuple):
    """:ChannelProfile command (RVP)."""

    name: str = "chn_XXX"
    bed_slope: float = 0
    survey_points: List[Tuple[float, float]] = []
    roughness_zones: List[Tuple[float, float]] = []

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
        d = self._asdict()
        d["survey_points"] = "\n".join(f"\t\t{p[0]} {p[1]}" for p in d["survey_points"])
        d["roughness_zones"] = "\n".join(
            f"\t\t{z[0]} {z[1]}" for z in d["roughness_zones"]
        )
        return pat.format(**d)


class GridWeightsCommand(NamedTuple):
    """:GridWeights command."""

    number_hrus: int = 0
    number_grid_cells: int = 0
    data: List[Tuple[int, int, float]] = []

    def to_rv(self, indent_level=0):
        indent = "\t" * indent_level
        pat = """
{indent}:GridWeights
{indent}\t:NumberHRUs {number_hrus}
{indent}\t:NumberGridCells {number_grid_cells}
{data}
{indent}:EndGridWeights
        """.strip()
        d = self._asdict()
        d["indent"] = indent
        d["data"] = "\n".join(f"{indent}\t{p[0]} {p[1]} {p[2]}" for p in self.data)
        return pat.format(**d)


class GriddedForcingCommand(NamedTuple):
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
\t:DimNamesNC {dim_names_nc}
{grid_weights}
:EndGriddedForcing
        """
        d = self._asdict()
        d["dim_names_nc"] = " ".join(self.dim_names_nc)
        d["grid_weights"] = self.grid_weights.to_rv()
        return pat.format(**d)
