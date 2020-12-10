from dataclasses import dataclass
from typing import List, NamedTuple, Tuple


# @dataclass
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


# @dataclass
class SubBasinsCommand(NamedTuple):
    """RVH :SubBasins command."""

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
    """Record to populate RVH :HRUs command internal table."""

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
    """RVH :HRUs command."""

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
    """RVH :Reservoir command."""

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
    """RVH :SubBasinGroup command."""

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
    """RVP :ChannelProfile command."""

    name: str = "chn_XXX"
    bed_slope: float = 0
    survey_points: List[Tuple[float, float]] = []
    roughness_zones: List[Tuple[float, float]] = []

    def to_rv(self):
        pat = """
:ChannelProfile	{name}
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
