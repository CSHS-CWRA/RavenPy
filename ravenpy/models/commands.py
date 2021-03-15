import re
from dataclasses import asdict, dataclass, field
from textwrap import dedent
from typing import Any, Dict, Optional, Tuple

INDENT = " " * 4
VALUE_PADDING = 10


class RavenConfig:
    def __str__(self):
        return self.to_rv()


@dataclass
class SubBasinsCommand(RavenConfig):
    """SubBasins command (RVH)."""

    @dataclass
    class Record(RavenConfig):
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
            return " ".join(f"{v: <{VALUE_PADDING}}" for v in d.values())

    subbasins: Tuple[Record] = ()

    template = """
    :SubBasins
        :Attributes   ID NAME DOWNSTREAM_ID PROFILE REACH_LENGTH  GAUGED
        :Units      none none          none    none           km    none
    {subbasin_records}
    :EndSubBasins
    """

    def to_rv(self):
        recs = [f"    {sb}" for sb in self.subbasins]
        return dedent(self.template).format(subbasin_records="\n".join(recs))


@dataclass
class HRUsCommand(RavenConfig):
    """HRUs command (RVH)."""

    @dataclass
    class Record(RavenConfig):
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

        def to_rv(self):
            d = asdict(self)
            return " ".join(f"{v: <{VALUE_PADDING * 2}}" for v in d.values())

    hrus: Tuple[Record] = ()

    template = """
    :HRUs
        :Attributes      AREA  ELEVATION       LATITUDE      LONGITUDE BASIN_ID       LAND_USE_CLASS           VEG_CLASS      SOIL_PROFILE  AQUIFER_PROFILE TERRAIN_CLASS      SLOPE     ASPECT
        :Units            km2          m            deg            deg     none                  none               none              none             none          none        deg       degN
    {hru_records}
    :EndHRUs
    """

    def to_rv(self):
        recs = [f"    {hru}" for hru in self.hrus]
        return dedent(self.template).format(hru_records="\n".join(recs))


@dataclass
class ReservoirCommand(RavenConfig):
    """Reservoir command (RVH)."""

    subbasin_id: int = 0
    hru_id: int = 0
    name: str = "Lake_XXX"
    weir_coefficient: float = 0
    crest_width: float = 0
    max_depth: float = 0
    lake_area: float = 0  # in m^2

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

    def to_rv(self):
        d = asdict(self)
        return dedent(self.template).format(**d)


@dataclass
class SubBasinGroupCommand(RavenConfig):
    """SubBasinGroup command (RVH)."""

    name: str = ""
    subbasin_ids: Tuple[int] = ()

    template = """
    :SubBasinGroup {name}
        {subbasin_ids}
    :EndSubBasinGroup
    """

    def to_rv(self):
        d = asdict(self)
        n_per_line = 10
        sbids = sorted(self.subbasin_ids)
        sbids_lines = [
            map(str, sbids[i : i + n_per_line])
            for i in range(0, len(sbids), n_per_line)
        ]
        d["subbasin_ids"] = "\n    ".join([" ".join(sbids) for sbids in sbids_lines])
        return dedent(self.template).format(**d)


@dataclass
class SBGroupPropertyMultiplierCommand(RavenConfig):

    group_name: str = ""
    parameter_name: str = ""
    mult: float = None

    template = ":SBGroupPropertyMultiplier {group_name} {parameter_name} {mult}"

    def to_rv(self):
        return dedent(self.template).format(**asdict(self))


@dataclass
class ChannelProfileCommand(RavenConfig):
    """ChannelProfile command (RVP)."""

    name: str = "chn_XXX"
    bed_slope: float = 0
    survey_points: Tuple[Tuple[float, float]] = ()
    roughness_zones: Tuple[Tuple[float, float]] = ()

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

    def to_rv(self):
        d = asdict(self)
        d["survey_points"] = "\n".join(
            f"{INDENT * 2}{p[0]} {p[1]}" for p in d["survey_points"]
        )
        d["roughness_zones"] = "\n".join(
            f"{INDENT * 2}{z[0]} {z[1]}" for z in d["roughness_zones"]
        )
        return dedent(self.template).format(**d)


@dataclass
class GridWeightsCommand(RavenConfig):
    """GridWeights command."""

    number_hrus: int = 0
    number_grid_cells: int = 0
    data: Tuple[Tuple[int, int, float]] = ()

    template = """
    {indent}:GridWeights
    {indent}    :NumberHRUs {number_hrus}
    {indent}    :NumberGridCells {number_grid_cells}
    {data}
    {indent}:EndGridWeights
    """

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
        n_hrus, n_grid_cells, data = m.groups()
        data = [d.strip().split() for d in data.split("\n")]
        data = tuple((int(h), int(c), float(w)) for h, c, w in data)
        return cls(
            number_hrus=int(n_hrus), number_grid_cells=int(n_grid_cells), data=data
        )

    def to_rv(self, indent_level=0):
        indent = INDENT * indent_level
        d = asdict(self)
        d["indent"] = indent
        d["data"] = "\n".join(f"{indent}    {p[0]} {p[1]} {p[2]}" for p in self.data)
        return dedent(self.template).strip().format(**d)


@dataclass
class GriddedForcingCommand(RavenConfig):
    """GriddedForcing command (RVT)."""

    name: str = ""
    forcing_type: str = ""
    file_name_nc: str = ""
    var_name_nc: str = ""
    dim_names_nc: Tuple[str, str, str] = ("x", "y", "t")
    time_shift: Optional[int] = None  # in days
    linear_transform: Optional[Tuple[float, float]] = None
    grid_weights: GridWeightsCommand = None

    template = """
    :GriddedForcing {name}
        :ForcingType {forcing_type}
        :FileNameNC {file_name_nc}
        :VarNameNC {var_name_nc}
        :DimNamesNC {dim_names_nc}
        {time_shift}
        {linear_transform}
        {grid_weights}
    :EndGriddedForcing
    """

    def to_rv(self):
        d = asdict(self)
        d["dim_names_nc"] = " ".join(self.dim_names_nc)
        d["time_shift"] = f":TimeShift {self.time_shift}" if self.time_shift else ""
        if self.linear_transform:
            slope, intercept = self.linear_transform
            d["linear_transform"] = f":LinearTransform {slope} {intercept}"
        else:
            d["linear_transform"] = ""
        d["grid_weights"] = self.grid_weights.to_rv(indent_level=1)
        return dedent(self.template).format(**d)


@dataclass
class BaseValueCommand(RavenConfig):
    """BaseValueCommand."""

    tag: str = ""
    value: Any = None

    template = ":{tag} {value}"

    # Overloading init to freeze the tag.
    def __init__(self, value):
        self.value = value

    def to_rv(self):
        return self.template.format(**asdict(self))


@dataclass
class BaseBooleanCommand(RavenConfig):
    tag: str = ""
    value: bool = False

    template = ":{tag}"

    def to_rv(self):
        return self.template.format(tag=self.tag) if self.value else ""


@dataclass
class LinearTransform(RavenConfig):
    scale: float = 1
    offset: float = 0

    template = ":LinearTransform {slope:.15f} {intercept:.15f}"

    def to_rv(self):
        return self.template.format(**asdict(self))


@dataclass
class DataCommand(RavenConfig):
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

    template = """
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
    :End{kind}
    """

    def to_rv(self):
        d = asdict(self)
        d["kind"] = (
            "ObservationData"
            if self.var == "water_volume_transport_in_river_channel"
            else "Data"
        )
        return dedent(self.template).format(**d)


@dataclass
class HRUStateVariableTableCommand(RavenConfig):
    """Initial condition for a given HRU."""

    @dataclass
    class Record(RavenConfig):
        index: int = 1
        surface_water: float = 0
        atmosphere: float = 0
        atmos_precip: float = 0
        ponded_water: float = 0
        soil0: float = 0
        soil1: float = 0
        soil2: float = 0
        soil3: float = 0
        snow_temp: float = 0
        snow: float = 0
        snow_cover: float = 0
        aet: float = 0
        convolution0: float = 0
        convolution1: float = 0
        conv_stor0: float = 0
        conv_stor1: float = 0
        conv_stor2: float = 0
        conv_stor3: float = 0
        conv_stor4: float = 0
        conv_stor5: float = 0
        conv_stor6: float = 0
        conv_stor7: float = 0
        conv_stor8: float = 0
        conv_stor9: float = 0
        conv_stor10: float = 0
        conv_stor11: float = 0
        conv_stor12: float = 0
        conv_stor13: float = 0
        conv_stor14: float = 0
        conv_stor15: float = 0
        conv_stor16: float = 0
        conv_stor17: float = 0
        conv_stor18: float = 0
        conv_stor19: float = 0
        conv_stor20: float = 0
        conv_stor21: float = 0
        conv_stor22: float = 0
        conv_stor23: float = 0
        conv_stor24: float = 0
        conv_stor25: float = 0
        conv_stor26: float = 0
        conv_stor27: float = 0
        conv_stor28: float = 0
        conv_stor29: float = 0
        conv_stor30: float = 0
        conv_stor31: float = 0
        conv_stor32: float = 0
        conv_stor33: float = 0
        conv_stor34: float = 0
        conv_stor35: float = 0
        conv_stor36: float = 0
        conv_stor37: float = 0
        conv_stor38: float = 0
        conv_stor39: float = 0
        conv_stor40: float = 0
        conv_stor41: float = 0
        conv_stor42: float = 0
        conv_stor43: float = 0
        conv_stor44: float = 0
        conv_stor45: float = 0
        conv_stor46: float = 0
        conv_stor47: float = 0
        conv_stor48: float = 0
        conv_stor49: float = 0
        conv_stor50: float = 0
        conv_stor51: float = 0
        conv_stor52: float = 0
        conv_stor53: float = 0
        conv_stor54: float = 0
        conv_stor55: float = 0
        conv_stor56: float = 0
        conv_stor57: float = 0
        conv_stor58: float = 0
        conv_stor59: float = 0
        conv_stor60: float = 0
        conv_stor61: float = 0
        conv_stor62: float = 0
        conv_stor63: float = 0
        conv_stor64: float = 0
        conv_stor65: float = 0
        conv_stor66: float = 0
        conv_stor67: float = 0
        conv_stor68: float = 0
        conv_stor69: float = 0
        conv_stor70: float = 0
        conv_stor71: float = 0
        conv_stor72: float = 0
        conv_stor73: float = 0
        conv_stor74: float = 0
        conv_stor75: float = 0
        conv_stor76: float = 0
        conv_stor77: float = 0
        conv_stor78: float = 0
        conv_stor79: float = 0
        conv_stor80: float = 0
        conv_stor81: float = 0
        conv_stor82: float = 0
        conv_stor83: float = 0
        conv_stor84: float = 0
        conv_stor85: float = 0
        conv_stor86: float = 0
        conv_stor87: float = 0
        conv_stor88: float = 0
        conv_stor89: float = 0
        conv_stor90: float = 0
        conv_stor91: float = 0
        conv_stor92: float = 0
        conv_stor93: float = 0
        conv_stor94: float = 0
        conv_stor95: float = 0
        conv_stor96: float = 0
        conv_stor97: float = 0
        conv_stor98: float = 0
        conv_stor99: float = 0

        def to_rv(self):
            return ",".join(map(repr, asdict(self).values()))

    template = """
    :HRUStateVariableTable
    :Attributes,SURFACE_WATER,ATMOSPHERE,ATMOS_PRECIP,PONDED_WATER,SOIL[0],SOIL[1],SOIL[2],SOIL[3],SNOW_TEMP,SNOW,SNOW_COVER,AET,CONVOLUTION[0],CONVOLUTION[1],CONV_STOR[0],CONV_STOR[1],CONV_STOR[2],CONV_STOR[3],CONV_STOR[4],CONV_STOR[5],CONV_STOR[6],CONV_STOR[7],CONV_STOR[8],CONV_STOR[9],CONV_STOR[10],CONV_STOR[11],CONV_STOR[12],CONV_STOR[13],CONV_STOR[14],CONV_STOR[15],CONV_STOR[16],CONV_STOR[17],CONV_STOR[18],CONV_STOR[19],CONV_STOR[20],CONV_STOR[21],CONV_STOR[22],CONV_STOR[23],CONV_STOR[24],CONV_STOR[25],CONV_STOR[26],CONV_STOR[27],CONV_STOR[28],CONV_STOR[29],CONV_STOR[30],CONV_STOR[31],CONV_STOR[32],CONV_STOR[33],CONV_STOR[34],CONV_STOR[35],CONV_STOR[36],CONV_STOR[37],CONV_STOR[38],CONV_STOR[39],CONV_STOR[40],CONV_STOR[41],CONV_STOR[42],CONV_STOR[43],CONV_STOR[44],CONV_STOR[45],CONV_STOR[46],CONV_STOR[47],CONV_STOR[48],CONV_STOR[49],CONV_STOR[50],CONV_STOR[51],CONV_STOR[52],CONV_STOR[53],CONV_STOR[54],CONV_STOR[55],CONV_STOR[56],CONV_STOR[57],CONV_STOR[58],CONV_STOR[59],CONV_STOR[60],CONV_STOR[61],CONV_STOR[62],CONV_STOR[63],CONV_STOR[64],CONV_STOR[65],CONV_STOR[66],CONV_STOR[67],CONV_STOR[68],CONV_STOR[69],CONV_STOR[70],CONV_STOR[71],CONV_STOR[72],CONV_STOR[73],CONV_STOR[74],CONV_STOR[75],CONV_STOR[76],CONV_STOR[77],CONV_STOR[78],CONV_STOR[79],CONV_STOR[80],CONV_STOR[81],CONV_STOR[82],CONV_STOR[83],CONV_STOR[84],CONV_STOR[85],CONV_STOR[86],CONV_STOR[87],CONV_STOR[88],CONV_STOR[89],CONV_STOR[90],CONV_STOR[91],CONV_STOR[92],CONV_STOR[93],CONV_STOR[94],CONV_STOR[95],CONV_STOR[96],CONV_STOR[97],CONV_STOR[98],CONV_STOR[99]
        :Units,mm,mm,mm,mm,mm,mm,mm,mm,C,mm,0-1,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm,mm
        {hru_states}
    :EndHRUStateVariableTable
    """
    hru_states: Dict[int, Record] = field(default_factory=dict)

    def to_rv(self):
        return dedent(self.template).format(
            hru_states="\n".join(map(str, self.hru_states.values()))
        )


@dataclass
class BasinIndexCommand(RavenConfig):
    """Initial conditions for a flow segment."""

    index: int = 1
    name: str = "watershed"
    channelstorage: float = 0
    rivuletstorage: float = 0
    qout: tuple = (0,)
    qoutlast: float = 0
    qlat: tuple = (0, 0, 0)
    qlatlast: float = 0
    qin: tuple = 20 * (0,)

    template = """
    :BasinIndex {index},{name}
        :ChannelStorage, {channelstorage}
        :RivuletStorage, {rivuletstorage}
        :Qout,{nsegs},{qout},{qoutlast}
        :Qlat,{nQlatHist},{qlat},{qlatlast}
        :Qin ,{nQinHist}, {qin}
        """

    def to_rv(self):
        return (
            dedent(self.template)
            .format(
                **asdict(self),
                nsegs=len(self.qout),
                nQlatHist=len(self.qlat),
                nQinHist=len(self.qin),
            )
            .replace("(", "")
            .replace(")", "")
        )


@dataclass
class BasinStateVariablesCommand(RavenConfig):

    basin_states: Dict[int, BasinIndexCommand] = field(default_factory=dict)

    template = """
    :BasinStateVariables
        {basin_states_list}
    :EndBasinStateVariables
    """

    def to_rv(self):
        return dedent(self.template).format(
            basin_states_list="\n".join(map(str, self.basin_states.values()))
        )


@dataclass
class SoilClassesCommand(RavenConfig):
    @dataclass
    class Record(RavenConfig):
        name: str = ""

        def to_rv(self):
            return " ".join(map(str, asdict(self).values()))

    soil_classes: Tuple[Record] = ()

    template = """
    :SoilClasses
        {soil_class_records}
    :EndSoilClasses
    """

    def to_rv(self):
        return dedent(self.template).format(
            soil_class_records="\n".join(map(str, self.soil_classes))
        )


@dataclass
class SoilProfilesCommand(RavenConfig):
    @dataclass
    class Record(RavenConfig):
        name: str = ""
        number_of_layers: int = 1
        soil_class: str = ""
        thickness: float = 0

        def to_rv(self):
            return " ".join(map(str, asdict(self).values()))

    soil_profiles: Tuple[Record] = ()

    template = """
    :SoilProfiles
        # name, number of layers, soil class, thickness [m]
        {soil_profile_records}
    :EndSoilProfiles
    """

    def to_rv(self):
        return dedent(self.template).format(
            soil_profile_records="\n".join(map(str, self.soil_profiles))
        )


@dataclass
class VegetationClassesCommand(RavenConfig):
    @dataclass
    class Record(RavenConfig):
        name: str = ""
        max_ht: float = 0
        max_lai: float = 0
        max_leaf_cond: float = 0

        def to_rv(self):
            return " ".join(map(str, asdict(self).values()))

    vegetation_classes: Tuple[Record] = ()

    template = """
    :VegetationClasses
        :Attributes,                MAX_HT,       MAX_LAI,    MAX_LEAF_COND
        :Units,                       m,            none,       mm_per_s
        {vegetation_class_records}
    :EndVegetationClasses
    """

    def to_rv(self):
        return dedent(self.template).format(
            vegetation_class_records="\n".join(map(str, self.vegetation_classes))
        )


@dataclass
class LandUseClassesCommand(RavenConfig):
    @dataclass
    class Record(RavenConfig):
        name: str = ""
        impermeable_frac: float = 0
        forest_coverage: float = 0

        def to_rv(self):
            return " ".join(map(str, asdict(self).values()))

    land_use_classes: Tuple[Record] = ()

    template = """
    :LandUseClasses
        :Attributes,        IMPERMEABLE_FRAC,         FOREST_COVERAGE
        :Units,                     fract,                    fract
        {land_use_class_records}
    :EndLandUseClasses
    """

    def to_rv(self):
        return dedent(self.template).format(
            land_use_class_records="\n".join(map(str, self.land_use_classes))
        )


@dataclass
class ObservationDataCommand(RavenConfig):

    data_type: str = "HYDROGRAPH"
    subbasin_id: int = None
    units: str = "m3/s"
    file_name_nc: str = None
    var_name_nc: str = None
    dim_names_nc: Tuple[str, str] = ("?", "?")
    station_idx: int = None

    template = """
    :ObservationData {data_type} {subbasin_id} {units}
        :ReadFromNetCDF
             :FileNameNC     {file_name_nc}
             :VarNameNC      {var_name_nc}
             :DimNamesNC     {dim_names_nc}
             :StationIdx     {station_idx}
        :EndReadFromNetCDF
    :EndObservationData
    """

    def to_rv(self):
        d = asdict(self)
        dns = d["dim_names_nc"]
        d["dim_names_nc"] = f"{dns[0]} {dns[1]}"
        return dedent(self.template).format(**d)


class RainCorrection(BaseValueCommand):
    tag: str = "RainCorrection"
    value: float = 1.0


class SnowCorrection(BaseValueCommand):
    tag: str = "SnowCorrection"
    value: float = 1.0
