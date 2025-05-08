import datetime as dt
import itertools
import logging
import re
from collections.abc import Sequence
from pathlib import Path
from textwrap import dedent, indent
from typing import (
    Annotated,
    Any,
    Literal,
    Optional,
    Union,
    get_args,
    get_origin,
    no_type_check,
)

import cftime  # noqa: F401
import xarray as xr
from pydantic import (
    ConfigDict,
    Field,
    FilePath,
    HttpUrl,
    NonNegativeFloat,
    PositiveInt,
    PrivateAttr,
    ValidationError,
    field_validator,
    model_validator,
)

from ..config import options
from .base import (
    Command,
    FlatCommand,
    GenericParameterList,
    LineCommand,
    ListCommand,
    ParameterList,
    Record,
    RootCommand,
    RootRecord,
    Sym,
)
from .utils import filter_for, get_annotations, nc_specs

LOGGER = logging.getLogger("RavenPy")


"""
Raven commands
==============

The syntax of those commands match as closely as possible the Raven documentation.
"""


# TODO: Strip command outputs of white space to facilitate testing.
# TODO: Add docstrings
# TODO: Create tests in tests/config/test_commands.py


INDENT = " " * 4
VALUE_PADDING = 10
T12 = tuple[
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


class Process(Command):
    """
    Process type embedded in HydrologicProcesses command.

    See processes.py for list of processes.
    """

    algo: str = "RAVEN_DEFAULT"
    source: Optional[str] = None
    to: Sequence[str] = ()

    _sub = ""

    @field_validator("to", mode="before")
    @classmethod
    def is_list(cls, v):
        if isinstance(v, str):
            return [
                v,
            ]
        return v

    @field_validator("source")
    @classmethod
    def is_source_state_variable(cls, v: str):
        if v is not None:
            if re.split(r"\[", v)[0] not in get_args(options.StateVariables):
                raise ValueError(f"{v} not recognized as a state variable.")
        return v

    @field_validator("to")
    @classmethod
    def is_to_state_variable(cls, v: Sequence):
        for x in v:
            if re.split(r"\[", x)[0] not in get_args(options.StateVariables):
                raise ValueError(f"{x} not recognized as a state variable.")
        return v

    def to_rv(self):
        t = " ".join(map("{:<19}".format, self.to))
        p = getattr(self, "p", None)
        p = p if p is not None else ""
        cmd = self._sub + self.__class__.__name__
        return f":{cmd:<20} {self.algo:<19} {self.source:<19} {t} {p}".strip()


class SoilModel(RootCommand):
    """:SoilModel SOIL_MULTILAYER 6"""

    root: int

    def to_rv(self):
        cmd = "SoilModel"
        return f":{cmd:<20} SOIL_MULTILAYER {self.root}\n"


class LinearTransform(Command):
    """:LinearTransform 1.0 -273.15"""

    scale: float = 1
    offset: float = 0

    def to_rv(self):
        cmd = "LinearTransform"
        if (self.scale != 1) or (self.offset != 0):
            return f":{cmd:<20} {self.scale:.15f} {self.offset:.15f}\n"
        return ""


class RainSnowTransition(Command):
    """Specify the range of temperatures over which there will be a rain/snow mix when partitioning total precipitation into rain/snow components.

    :RainSnowTransition [temp] [delta]

    # equivalent to (the preferred option)
    :GlobalParameter RAINSNOW_TEMP [rainsnow_temp]
    :GlobalParameter RAINSNOW_DELTA [rainsnow_delta]
    """

    temp: Sym
    """Midpoint of the temperature range [C]."""
    delta: Sym
    """Range [C]."""

    def to_rv(self):
        cmd = "RainSnowTransition"
        return f":{cmd:<20} {self.temp} {self.delta}\n"


class EvaluationPeriod(LineCommand):
    """:EvaluationPeriod [period_name] [start yyyy-mm-dd] [end yyyy-mm-dd]"""

    name: str
    start: dt.date
    end: dt.date


class CustomOutput(LineCommand):
    """Create custom output file to track a single variable, parameter, or forcing function over time at a number of basins, HRUs, or across watershed.

    :CustomOutput DAILY AVERAGE AET BY_HRU
    """  # noqa: E501

    time_per: Literal["DAILY", "MONTHLY", "YEARLY", "WATER_YEARLY", "CONTINUOUS"]
    """Time period."""
    stat: Literal["AVERAGE", "MAXIMUM", "MINIMUM", "RANGE", "MEDIAN", "QUARTILES"]
    """Statistic reported for each time interval."""
    variable: str
    """
    Variable or parameter name.

    Consult the Raven documentation for the list of allowed names.
    """
    space_agg: Literal[
        "BY_BASIN", "BY_HRU", "BY_HRU_GROUP", "BY_SB_GROUP", "ENTIRE_WATERSHED"
    ]
    """Spatial evaluation domain."""
    filename: str = ""
    """
    Output file name.

    Defaults to something approximately like: `<run name>_<variable>_<time_per>_<stat>_<space_agg>.nc`.
    """


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
    """SoilProfiles command.

    Example::

        :SoilProfiles
          # name, #horizons, hor1, th1, hor2, th2
          LAKE, 0
          GLACIER, 0
          LOAM_SEQ, 2, LOAM, 0.5, SAND, 1.5
          ALL_SAND, 2, SAND, 0.5, SAND, 1.5
        :EndSoilProfiles
    """

    root: Sequence[SoilProfile]


class SubBasin(Record):
    """Record to populate RVH :SubBasins command internal table."""

    subbasin_id: PositiveInt = 1
    name: str = "sub_001"
    downstream_id: int = -1
    profile: str = "NONE"
    reach_length: float | str = 0
    gauged: bool = True
    gauge_id: Optional[str] = ""  # This attribute is not rendered to RVH

    @classmethod
    @field_validator("reach_length")
    def check_reach_length(cls, v):
        if v == "ZERO-":
            return 0
        return v

    def __str__(self):
        d = self.model_dump()
        d["gauged"] = int(d["gauged"])
        del d["gauge_id"]
        return " ".join(f"{v: <{VALUE_PADDING}}" for v in d.values())

    model_config = ConfigDict(coerce_numbers_to_str=True)


class SubBasins(ListCommand):
    """SubBasins command (RVH)

    Example::

        :SubBasins
          :Attributes, NAME, DOWNSTREAM_ID, PROFILE, REACH_LENGTH, GAUGED
          :Units, none, none, none, km, none
          1, Downstream, -1, DEFAULT, 3.0, 1
          2, Upstream, 1, DEFAULT, 3.0, 0
        :EndSubBasins
    """

    root: Sequence[SubBasin]
    _Attributes: Sequence[str] = PrivateAttr(
        ["ID", "NAME", "DOWNSTREAM_ID", "PROFILE", "REACH_LENGTH", "GAUGED"]
    )
    _Units: Sequence[str] = PrivateAttr(["none", "none", "none", "none", "km", "none"])

    @classmethod
    def parse(cls, s):
        """Parse a SubBasins command."""
        pat = r":SubBasins(.+):EndSubBasins" ""
        out = []
        keys = [
            "subbasin_id",
            "name",
            "downstream_id",
            "profile",
            "reach_length",
            "gauged",
        ]

        # Remove all commented lines (starting with #)
        s = re.sub(r"^\s*#.*$", "", s, flags=re.MULTILINE)
        if match := re.search(pat, s.strip(), re.DOTALL):
            spat = r"^\s*([^:^\s].+)$"
            for line in re.findall(spat, match.groups()[0], re.MULTILINE):
                # Split line into fields
                values = re.split(r"\s+", line.strip())
                # Convert to SubBasin record
                out.append(SubBasin(**dict(zip(keys, values))))
        return cls(root=out)


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
    aspect: Annotated[float, Field(ge=0, le=360)] = 0.0
    # This field is not part of the Raven config, it is needed for serialization,
    # to specify which HRU subclass to use when necessary
    hru_type: Optional[str] = None

    def __str__(self):
        import numpy as np

        d = self.model_dump()
        del d["hru_type"]

        # Adjust horizontal spacing to match attributes names
        attrs = ["HRU_ID"] + HRUs._Attributes.default
        ns = np.array(list(map(len, attrs)))
        ns[1:] += 1

        # In Py3.10 zip(strict=True) can be used instead.
        if len(d) != len(attrs):
            raise ValueError(
                f"HRU record has {len(d)} attributes, but {len(attrs)} are expected."
            )
        return " ".join(f"{v:<{n}}" for v, n in zip(d.values(), ns))


class HRUs(ListCommand):
    """HRUs command (RVH)."""

    root: Sequence[HRU]
    _Attributes: Sequence[str] = PrivateAttr(
        [
            "AREA",
            "ELEVATION",
            "LATITUDE",
            "LONGITUDE",
            "BASIN_ID",
            "LAND_USE_CLASS",
            "VEG_CLASS",
            "SOIL_PROFILE",
            "AQUIFER_PROFILE",
            "TERRAIN_CLASS",
            "SLOPE",
            "ASPECT",
        ]
    )
    _Units: Sequence[str] = PrivateAttr(
        [
            "km2",
            "m",
            "deg",
            "deg",
            "none",
            "none",
            "none",
            "none",
            "none",
            "none",
            "deg",
            "degN",
        ]
    )

    @field_validator("root", mode="before")
    @classmethod
    def ignore_unrecognized_hrus(cls, values):
        """Ignore HRUs with unrecognized hru_type.

        HRUs are ignored only if all allowed HRU classes define `hru_type`, and if the values passed include it.
        """
        import collections
        import warnings

        a = cls.model_fields["root"].annotation

        # Annotation should be a sequence
        if get_origin(a) != collections.abc.Sequence:
            return values

        # Extract allowed HRU types
        allowed = [hru.model_fields["hru_type"].default for hru in get_annotations(a)]

        # If some HRU classes do not define rhu_type, skip filtering
        if None in allowed:
            return values

        allowed.append(None)

        out = [value for value in values if getattr(value, "hru_type", None) in allowed]
        if len(out) != len(values):
            warnings.warn(
                "HRUs with an unrecognized `hru_type` attribute were ignored."
            )
        return out

    @classmethod
    def parse(cls, s):
        """Parse a HRUs command."""
        pat = r":HRUs(.+):EndHRUs"
        out = []
        keys = [
            "hru_id",
            "area",
            "elevation",
            "latitude",
            "longitude",
            "subbasin_id",
            "land_use_class",
            "veg_class",
            "soil_profile",
            "aquifer_profile",
            "terrain_class",
            "slope",
            "aspect",
        ]
        # Remove all commented lines (starting with #)
        s = re.sub(r"^\s*#.*$", "", s, flags=re.MULTILINE)

        if match := re.search(pat, s.strip(), re.DOTALL):
            spat = r"^\s*([^:^\s].+)$"
            for line in re.findall(spat, match.groups()[0], re.MULTILINE):
                # Split line into fields
                values = re.split(r"\s+", line.strip())
                # Convert to HRU record
                out.append(dict(zip(keys, values)))
        return cls(root=out)

    def rename(self, mapping):
        """Rename HRU attributes.

        Parameters
        ----------
        mapping : dict
            Nested dictionary keyed by HRU attributes, with keys as old names and values as new names.
        """
        out = self.copy()
        for attr, amap in mapping.items():
            for hru in out.root:
                cur = getattr(hru, attr)
                for old, new in amap.items():
                    if cur == old:
                        setattr(hru, attr, new)
        return out


class HRUGroup(FlatCommand):
    class _Rec(RootRecord):
        root: Sequence[str]

        def __str__(self):
            return ",".join(self.root)

    name: str
    groups: _Rec

    @property
    def _template(self):
        return """
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
    max_depth: float = Field(0, alias="MaxDepth", description="Max depth (m)")
    lake_area: float = Field(0, alias="LakeArea", description="Lake area (m2)")
    absolute_crest_height: Optional[float] = Field(
        None, alias="AbsoluteCrestHeight", description="Absolute crest height (m)"
    )
    max_capacity: Optional[float] = Field(
        None, alias="MaxCapacity", description="Maximum capacity in m3"
    )

    class SeepageParameters(LineCommand):
        """:SeepageParameters [K_seep] [href]"""

        k_seep: float
        h_ref: float

    seepage_parameters: Optional[SeepageParameters] = Field(
        None, alias="SeepageParameters"
    )

    class StageRelations(RootCommand):
        """Stage relations for the reservoir."""

        class StageRelation(RootRecord):
            """Stage relation record.

            h, q, v, a, [u]
            """

            root: (
                tuple[float, float, float, float]
                | tuple[float, float, float, float, float]
            ) = None

            def __str__(self):
                return ", ".join([f"{x:>12}" for x in self.root])

        root: Sequence[StageRelation] = Field(None)

    stage_relations: Optional[StageRelations] = Field(None, alias="StageRelations")

    class OutflowControlStructure(Command):
        """Outflow control structure for the reservoir."""

        target_subbasin_id: Optional[int] = Field(None, alias="TargetSubBasin")
        downstream_reference_elevation: Optional[float] = Field(
            None, alias="DownstreamReferenceElevation"
        )

        class StageDischargeTable(RootCommand):
            """Stage discharge table for the outflow control structure.

            Example::

                :StageDischargeTable C1 #one gate open
                      N
                      [h,Q]xN
                :EndStageDischargeTable
            """

            _n: Optional[int] = None

            class StageDischargeRecord(RootRecord):
                """Stage discharge record."""

                root: tuple[float, float] = None

            root: Sequence[StageDischargeRecord] = Field(None)

            @model_validator(mode="after")
            def _set_n(self):
                """Set the number of records in the table."""
                if self.root is not None and self.n is None:
                    self.n = len(self.root)
                return self

            @property
            def _template(self):
                return """
                :StageDischargeTable
                  {n}
                {_records}
                :EndStageDischargeTable
                """

        stage_discharge_table: Optional[Sequence[StageDischargeTable]] = Field(
            None, alias="StageDischargeTable"
        )

        class BasicWeir(LineCommand):
            """Basic weir for the outflow control structure.

            :BasicWeir [curve_name3] [elev] [crestwidth] [coeff]
            """

            curve_name: str
            elev: float
            crest_width: float
            coeff: float

        basic_weir: BasicWeir = Field(None, alias="BasicWeir")

        class OperatingRegime(Command):
            """Operating regime for the outflow control structure.

            :OperatingRegime [curve_name] [elev] [coeff]
            """

            name: str
            use_curve: str = Field(None, alias="UseCurve")
            condition: str = Field(None, alias="Condition")
            constraint: str = Field(None, alias="Constraint")

        operating_regime: OperatingRegime = Field(None, alias="OperatingRegime")

    @property
    def _template(self):
        return """
            :Reservoir {name}
            {_commands}
            :EndReservoir
            """

    @classmethod
    def parse(cls, s) -> list:
        """Return a list of Reservoir records parsed from an RV file."""
        # Remove all commented lines (starting with #)
        s = re.sub(r"^\s*#.*$", "", s, flags=re.MULTILINE)
        pat = r":Reservoir\s+(\w+)\s+(.+?):EndReservoir"
        out = []

        for name, content in re.findall(pat, s.strip(), re.DOTALL):
            # Extract parameters
            subbasin_id = re.search(r":SubBasinID\s+(\d+)", content).group(1)
            hru_id = re.search(r":HRUID\s+(\d+)", content).group(1)
            type_ = re.search(r":Type\s+(\w+)", content).group(1)
            weir_coefficient = re.search(
                r":WeirCoefficient\s+([\d.-]+)", content
            ).group(1)
            crest_width = re.search(r":CrestWidth\s+([\d.-]+)", content).group(1)
            max_depth = re.search(r":MaxDepth\s+([\d.-]+)", content).group(1)
            lake_area = re.search(r":LakeArea\s+([\d.-]+)", content).group(1)
            seepage_parameters = re.search(
                r":SeepageParameters\s+([\d.-]+)\s+([\d.-]+)", content
            ).groups()

            # Convert to Reservoir record
            out.append(
                cls(
                    name=name,
                    subbasin_id=subbasin_id,
                    hru_id=hru_id,
                    type=type_,
                    weir_coefficient=weir_coefficient,
                    crest_width=crest_width,
                    max_depth=max_depth,
                    lake_area=lake_area,
                    seepage_parameters=cls.SeepageParameters(
                        k_seep=seepage_parameters[0],
                        h_ref=seepage_parameters[1],
                    ),
                )
            )
        return out


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
        d = self.model_dump()
        n_per_line = 10
        sbids = sorted(self.sb_ids)
        sbids_lines = [
            map(str, sbids[i : i + n_per_line])
            for i in range(0, len(sbids), n_per_line)
        ]
        d["sb_ids"] = "\n  ".join([", ".join(sbids) for sbids in sbids_lines])
        return dedent(template).format(**d)

    @classmethod
    def parse(cls, s) -> list:
        """Parse a SubBasinGroup command and return a list of SubBasinGroup records."""
        # Remove all commented lines (starting with #)
        s = re.sub(r"^\s*#.*$", "", s, flags=re.MULTILINE)
        pat = r":SubBasinGroup\s+(\w+)\s+(.+?):EndSubBasinGroup"
        out = []

        for name, content in re.findall(pat, s.strip(), re.DOTALL):

            sb_ids = re.split(r"\s+", content.strip())

            # Convert to SubBasinGroup record
            out.append(SubBasinGroup(name=name, sb_ids=sb_ids))
        return out


class SBGroupPropertyMultiplier(LineCommand):
    """:SBGroupPropertyMultiplier [group_name] [parameter_name] [mult]"""

    group_name: str
    parameter_name: str
    mult: Sym


class ChannelProfile(FlatCommand):
    """ChannelProfile command (RVP)."""

    name: str = "chn_XXX"
    bed_slope: float = Field(0, alias="Bedslope")

    class SurveyPoints(ListCommand):
        """SurveyPoints

        [x, bed_elevation] x number of survey points.
        """

        class SurveyPoint(RootRecord):
            """SurveyPoint record."""

            root: tuple[float, float] = ()

        root: Sequence[SurveyPoint] = Field((SurveyPoint(),))

    survey_points: SurveyPoints = Field(SurveyPoints(), alias="SurveyPoints")

    class RoughnessZones(ListCommand):
        """RoughnessZones record.

        [x_zone, mannings_n] x number of roughness zones.
        """

        class RoughnessZone(RootRecord):
            """RoughnessZone record."""

            root: tuple[float, float] = ()

        root: Sequence[RoughnessZone] = Field((RoughnessZone(),))

    roughness_zones: RoughnessZones = Field(RoughnessZones(), alias="RoughnessZones")

    @property
    def _template(self):
        return """
               :{_cmd} {name}
               {_commands}{_records}
               :End{_cmd}
               """

    @classmethod
    def parse(cls, s) -> list:
        """Parse ChannelProfile commands and return a list of ChannelProfile records."""
        # Remove all commented lines (starting with #)
        s = re.sub(r"^\s*#.*$", "", s, flags=re.MULTILINE)

        pat = r":ChannelProfile\s+(\w+)\s+(.+?):EndChannelProfile"
        out = []

        for name, content in re.findall(pat, s.strip(), re.DOTALL):
            bed_slope = re.search(r":Bedslope\s+([\d.-]+)", content).group(1)
            survey_points = re.search(
                r":SurveyPoints(.+):EndSurveyPoints", content, re.DOTALL
            ).group(1)
            sp = [
                line.strip().split()
                for line in survey_points.splitlines()
                if line.strip()
            ]

            roughness_zones = re.search(
                r":RoughnessZones(.+):EndRoughnessZones", content, re.DOTALL
            ).group(1)
            rz = [line.split() for line in roughness_zones.splitlines() if line.strip()]
            # rz = re.split(r"\s+", roughness_zones.strip())

            # Convert to ChannelProfile record
            out.append(
                cls(
                    name=name,
                    bed_slope=bed_slope,
                    survey_points=sp,
                    roughness_zones=rz,
                )
            )
        return out


class GridWeights(Command):
    """GridWeights command.

    Notes
    -----
    command can be embedded in both a `GriddedForcing` or a `StationForcing`.

    The default is to have a single cell that covers an entire single HRU, with a
    weight of 1.
    """

    number_hrus: int = Field(1, alias="NumberHRUs")
    number_grid_cells: int = Field(1, alias="NumberGridCells")

    class GWRecord(RootRecord):
        """GridWeights record."""

        root: tuple[int, int, float] = (1, 0, 1.0)

    data: Sequence[GWRecord] = Field((GWRecord(),))

    # @field_validator("data")
    # @classmethod
    # def _sum_weights(cls, values):
    #     """Check that for each HRU weights sum to 1."""
    #     for hru, group in itertools.groupby(values, lambda x: x.root[0]):
    #         total = sum([x.root[2] for x in group])
    #         if not (0.999 < total < 1.001):
    #             raise ValueError(
    #                 f"GridWeights for HRU {hru} do not sum to 1.0: {total:.3f}"
    #             )
    #     return values

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
        n_hrus, n_grid_cells, data = m.groups()
        data = [d.strip().split() for d in data.split("\n")]
        data = tuple((int(h), int(c), float(w)) for h, c, w in data)
        return cls(
            number_hrus=int(n_hrus), number_grid_cells=int(n_grid_cells), data=data
        )


class RedirectToFile(RootCommand):
    """RedirectToFile command (RVT).

    Notes
    -----
    For the moment, this command can only be used in the context of a `GriddedForcingCommand`
    or a `StationForcingCommand`, as a `grid_weights` field replacement when inlining is not
    desired.
    """

    root: FilePath

    def to_rv(self):
        cmd = "RedirectToFile"
        return f":{cmd:<20} {self.root}\n"


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
    time_shift: Optional[float] = Field(
        None, alias="TimeShift", description="Time stamp shift in days."
    )
    linear_transform: Optional[LinearTransform] = Field(None, alias="LinearTransform")
    deaccumulate: Optional[bool] = Field(None, alias="Deaccumulate")

    latitude_var_name_nc: Optional[str] = Field(None, alias="LatitudeVarNameNC")
    longitude_var_name_nc: Optional[str] = Field(None, alias="LongitudeVarNameNC")
    elevation_var_name_nc: Optional[str] = Field(None, alias="ElevationVarNameNC")

    @field_validator("dim_names_nc")
    @classmethod
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

    @classmethod
    def from_nc(
        cls, fn, data_type, station_idx=None, alt_names=(), engine="h5netcdf", **kwds
    ):
        """Instantiate class from netCDF dataset."""
        specs = nc_specs(
            fn,
            data_type,
            station_idx=station_idx,
            alt_names=alt_names,
            engine=engine,
        )
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
    forcing_type: Optional[options.Forcings] = Field(None, alias="ForcingType")
    grid_weights: Union[GridWeights, RedirectToFile] = Field(
        GridWeights(), alias="GridWeights"
    )
    # StationIdx is not relevant to GriddedForcing
    station_idx: Optional[int] = None

    @property
    def _template(self):
        return """
        :{_cmd} {name}
        # HRU GridCell Weight
        {_commands}
        :End{_cmd}
        """

    @field_validator("dim_names_nc")
    @classmethod
    def check_dims(cls, v):
        if len(v) != 3:
            raise ValueError(
                "GriddedForcing netCDF datasets should have  three dimensions (lon, lat, time)."
            )
        return v


class StationForcing(GriddedForcing):
    """StationForcing command (RVT)."""

    @field_validator("dim_names_nc")
    @classmethod
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

    @property
    def _template(self):
        return """
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
    elevation: Optional[float] = Field(None, alias="Elevation")

    rain_correction: Optional[Sym] = Field(
        None, alias="RainCorrection", description="Rain correction"
    )
    snow_correction: Optional[Sym] = Field(
        None, alias="SnowCorrection", description="Snow correction"
    )

    monthly_ave_evaporation: Optional[Sequence] = Field(
        None, alias="MonthlyAveEvaporation"
    )
    monthly_ave_temperature: Optional[Sequence] = Field(
        None, alias="MonthlyAveTemperature"
    )
    monthly_min_temperature: Optional[Sequence] = Field(
        None, alias="MonthlyMinTemperature"
    )
    monthly_max_temperature: Optional[Sequence] = Field(
        None, alias="MonthlyMaxTemperature"
    )

    data: Optional[Sequence[Data]] = Field(None, alias="Data")

    _nested: bool = False

    @property
    def _template(self):
        return """
        :Gauge {name}
        {_commands}
        :EndGauge
        """

    @field_validator("monthly_ave_evaporation", "monthly_ave_temperature")
    @classmethod
    def confirm_monthly(cls, v):
        if v is not None and len(v) != 12:
            raise ValidationError("One value per month needed.")
        return v

    @classmethod
    def from_nc(
        cls,
        fn: Union[str, Path, Sequence[Path]],
        data_type: Optional[Sequence[str]] = None,
        station_idx: int = 1,
        alt_names: Optional[dict[str, str]] = None,
        mon_ave: bool = False,
        data_kwds: Optional[dict[str, Any]] = None,
        engine: str = "h5netcdf",
        **kwds,
    ) -> "Gauge":
        r"""Return Gauge instance with configuration options inferred from the netCDF itself.

        Parameters
        ----------
        fn : str or Path or Sequence[Path]
            NetCDF file path or paths.
        data_type : Sequence[str], optional
            Raven data types to extract from netCDF files, e.g. 'PRECIP', 'AVE_TEMP'. The algorithm tries to find all
            forcings in each file until one is found, then it stops searching for it in the following files.
        station_idx : int
            Index along station dimension. Starts at 1. Should be the same for all netCDF files.
        alt_names : dict
            Alternative variable names keyed by data type.
            Use this if variables do not correspond to CF standard defaults.
        mon_ave : bool
            If True, compute the monthly average.
        data_kwds : dict[options.Forcings, dict[str, str]]
            Additional `:Data` parameters keyed by forcing type and station id. Overrides inferred parameters.
            Use keyword "ALL" to pass parameters to all variables.
        engine : {"h5netcdf", "netcdf4", "pydap"}
            The engine used to open the dataset. Default is 'h5netcdf'.
        \*\*kwds : dict
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
                        fn=f,
                        data_type=dtype,
                        station_idx=idx,
                        alt_names=alt_names.get(dtype, ()),
                        mon_ave=mon_ave,
                        engine=engine,
                    )
                except ValueError as e:  # noqa: PERF203
                    LOGGER.info(e)
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
            raise ValueError(
                "No data found in netCDF files. Check that variable names follow CF conventions, "
                "or if not, provide `alt_names` mapping Raven data types to variable names."
            )

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


class ObservationData(Data, coerce_numbers_to_str=True):
    uid: str = "1"  # Basin ID or HRU ID
    data_type: Literal["HYDROGRAPH"] = "HYDROGRAPH"

    @property
    def _template(self):
        return """
        :ObservationData {data_type} {uid} {units}
        {_commands}
        :EndObservationData
        """

    @classmethod
    def from_nc(cls, fn, station_idx: int = 1, alt_names=(), engine="h5netcdf", **kwds):
        specs = nc_specs(
            fn,
            "HYDROGRAPH",
            station_idx=station_idx,
            alt_names=alt_names,
            engine=engine,
        )
        attrs = filter_for(cls, specs, **kwds, data_type="HYDROGRAPH")
        return cls(**attrs)


class HRUState(Record):
    hru_id: int = 1
    data: dict[str, Sym] = Field(default_factory=dict)

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

    root: Sequence[HRUState]
    _Attributes: Sequence[str] = PrivateAttr(None)

    @model_validator(mode="after")
    def set_attributes(self):
        from itertools import chain

        # Get all attribute names from HRUStates
        names = (
            sorted(list(set(chain(*[tuple(s.data.keys()) for s in self.root])))) or None
        )

        self.root = [
            HRUState(hru_id=s.hru_id, data={n: s.data.get(n, 0.0) for n in names})
            for s in self.root
        ]

        # Store names in private class attribute
        self._Attributes = names
        return self

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
        return cls(out)


class BasinIndex(Command, coerce_numbers_to_str=True):
    """Initial conditions for a flow segment."""

    # TODO: Check that qout and cie can be written with separating commas.
    sb_id: int = 1
    name: str = "watershed"
    channel_storage: float = Field(0.0, alias="ChannelStorage")
    rivulet_storage: float = Field(0.0, alias="RivuletStorage")
    qout: Sequence[float] = Field((1.0, 0.0, 0.0), alias="Qout")
    qlat: Optional[Sequence[float]] = Field(None, alias="Qlat")
    qin: Optional[Sequence[float]] = Field(None, alias="Qin")

    # model_config = ConfigDict(coerce_numbers_to_str=True)

    @property
    def _template(self):
        return """
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
                assert len(values) == 1  # noqa: S101
                rec_values[cmd[1:]] = float(values[0])
            else:
                rec_values[cmd[1:]] = tuple(values)
        return cls(**rec_values)


class BasinStateVariables(ListCommand):
    root: Sequence[BasinIndex]

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
        return cls(bs)


class SoilClasses(ListCommand):
    """SoilClasses command.

    Example::

        :SoilClasses
          :Attributes, %SAND, %CLAY, %SILT, %ORGANIC
          :Units, none, none, none, none
          SAND, 1, 0, 0, 0
          LOAM, 0.5, 0.1, 0.4, 0.4
        :EndSoilClasses
    """

    class SoilClass(Record):
        """SoilClass."""

        name: str
        mineral: Optional[tuple[float, float, float]] = None
        organic: Optional[float] = None

        @field_validator("mineral")
        @classmethod
        def validate_mineral(cls, v):
            """Assert sum of mineral fraction is 1."""
            if v is not None:
                # FIXME: assertions should not be found outside of testing code. Replace with conditional logic.
                assert sum(v) == 1, "Mineral fraction should sum to 1."  # noqa: S101
            return v

        @field_validator("mineral")
        @classmethod
        def validate_mineral_pct(cls, v: tuple):
            if v is not None:
                for x in v:
                    # FIXME: assertions should not be found outside of testing code. Replace with conditional logic.
                    assert (x >= 0) and (  # noqa: S101
                        x <= 1
                    ), "Value should be in [0,1]."
            return v

        @field_validator("organic")
        @classmethod
        def validate_organic_pct(cls, v: float):
            if v is not None:
                # FIXME: assertions should not be found outside of testing code. Replace with conditional logic.
                assert (v >= 0) and (v <= 1), "Value should be in [0,1]."  # noqa: S101
            return v

        def __str__(self):
            if self.mineral is None and self.organic is None:
                return self.name
            else:
                fmt = "{name:<16},{mineral[0]:>18},{mineral[1]:>18},{mineral[2]:>18},{organic:>18}"
                return fmt.format(**self.model_dump())

    root: Sequence[SoilClass] = ()
    _Attributes: Sequence[str] = PrivateAttr(["%SAND", "%CLAY", "%SILT", "%ORGANIC"])
    _Units: Sequence[str] = PrivateAttr(["none", "none", "none", "none"])


class VegetationClass(Record):
    name: str = ""
    max_ht: Sym = 0.0
    max_lai: Sym = 0.0
    max_leaf_cond: Sym = 0.0

    def __str__(self):
        template = "{name:<16},{max_ht:>18},{max_lai:>18},{max_leaf_cond:>18}"
        return template.format(**self.model_dump())


class VegetationClasses(ListCommand):
    root: Sequence[VegetationClass]
    _Attributes: Sequence[str] = PrivateAttr(["MAX_HT", "MAX_LAI", "MAX_LEAF_COND"])
    _Units: Sequence[str] = PrivateAttr(["m", "none", "mm_per_s"])


class LandUseClass(Record):
    name: str = ""
    impermeable_frac: Sym = 0.0
    forest_coverage: Sym = 0.0

    def __str__(self):
        template = "{name:<16},{impermeable_frac:>18},{forest_coverage:>18}"
        return template.format(**self.model_dump())


class LandUseClasses(ListCommand):
    root: Sequence[LandUseClass]
    _Attributes: Sequence[str] = PrivateAttr(["IMPERMEABLE_FRAC", "FOREST_COVERAGE"])
    _Units: Sequence[str] = PrivateAttr(["fract", "fract"])


class TerrainClass(Record):
    name: str
    hillslope_length: Sym
    drainage_density: Sym
    topmodel_lambda: Optional[Sym] = None

    def __str__(self):
        out = f"{self.name:<16},{self.hillslope_length:>18},{self.drainage_density:>18}"
        if self.topmodel_lambda is not None:
            out += f",{self.topmodel_lambda:>18}"
        return out


class TerrainClasses(ListCommand):
    root: Sequence[TerrainClass]
    _Attributes: Sequence[str] = PrivateAttr(
        ["HILLSLOPE_LENGTH", "DRAINAGE_DENSITY", "TOPMODEL_LAMBDA"]
    )
    _Units: Sequence[str] = PrivateAttr(["m", "km/km2"])


class _MonthlyRecord(Record):
    name: str = "[DEFAULT]"
    values: Sequence[float] = 12 * [1.0]

    def __str__(self):
        return f"{self.name:<16} " + ", ".join(map(str, self.values))


class SeasonalRelativeLAI(ListCommand):
    root: Sequence[_MonthlyRecord] = (_MonthlyRecord(),)


class SeasonalRelativeHeight(ListCommand):
    root: Sequence[_MonthlyRecord] = (_MonthlyRecord(),)


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
            "  ",
        )

    model_config = ConfigDict(coerce_numbers_to_str=True)


class SubBasinProperties(Command):
    parameters: Optional[Sequence[options.SubBasinProperties]] = Field(
        None, alias="Parameters"
    )
    records: Optional[Sequence[SubBasinProperty]] = Field(None)


class SoilParameterList(GenericParameterList):
    parameters: Optional[Sequence[options.SoilParameters]] = Field(
        None, alias="Parameters"
    )


class VegetationParameterList(GenericParameterList):
    parameters: Optional[Sequence[options.VegetationParameters]] = Field(
        None, alias="Parameters"
    )


class LandUseParameterList(GenericParameterList):
    parameters: Optional[Sequence[options.LandUseParameters]] = Field(
        None, alias="Parameters"
    )


class EnsembleMode(LineCommand):
    mode: Literal["ENSEMBLE_ENKF"] = "ENSEMBLE_ENKF"
    n: int = Field(..., description="Number of members")


class ObservationErrorModel(LineCommand):
    state: Literal["STREAMFLOW"]
    dist: Literal["DIST_UNIFORM", "DIST_NORMAL", "DIST_GAMMA"]
    p1: float
    p2: float
    adj: Literal["ADDITIVE", "MULTIPLICATIVE"]


class ForcingPerturbation(LineCommand):
    forcing: options.Forcings
    dist: Literal["DIST_UNIFORM", "DIST_NORMAL", "DIST_GAMMA"]
    p1: float
    p2: float
    adj: Literal["ADDITIVE", "MULTIPLICATIVE"]
    hru_grp: str = ""


class AssimilatedState(LineCommand):
    state: Union[options.StateVariables, Literal["STREAMFLOW"]]
    group: str


class AssimilateStreamflow(LineCommand, coerce_numbers_to_str=True):
    """Subbasin ID to assimilate streamflow for."""

    sb_id: str


# TODO: Harmonize this...
# Aliases for convenience
SC = SoilClasses.SoilClass
LU = LandUseClass
VC = VegetationClass
TC = TerrainClass
SB = SubBasin
SP = SoilProfile
PL = ParameterList
