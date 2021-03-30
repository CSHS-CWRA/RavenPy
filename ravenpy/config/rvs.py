from dataclasses import replace
from pathlib import Path
from textwrap import dedent

from ravenpy.config.commands import (
    DataCommand,
    GaugeCommand,
    GridWeightsCommand,
    HRUsCommand,
    SubBasinGroupCommand,
    SubBasinsCommand,
)


class RVH:
    def __init__(self, hrus, subbasins):
        self.hrus = hrus
        self.subbasins = subbasins
        self.land_subbasins: Tuple[int] = ()
        self.land_subbasin_property_multiplier = ""
        self.lake_subbasins = ()
        self.lake_subbasin_property_multiplier = ""
        self.reservoirs = ()

    def update(self, key, value):
        a = getattr(self, key, None)
        if a:
            setattr(self, key, value)
            return True
        return False

    def to_rv(self):
        tmpl = """
    {subbasins_cmd}

    {hrus_cmd}

    {land_subbasin_group_cmd}

    {lake_subbasin_group_cmd}

    {reservoir_cmd_list}
    """
        d = {
            "subbasins_cmd": SubBasinsCommand(self.subbasins),
            "hrus_cmd": HRUsCommand(self.hrus),
            "land_subbasin_group_cmd": SubBasinGroupCommand(
                "Land", self.land_subbasins
            ),
            "lake_subbasin_group_cmd": SubBasinGroupCommand(
                "Lakes", self.land_subbasins
            ),
            "reservoir_cmd_list": "\n\n".join(map(str, self.reservoirs)),
        }
        return dedent(tmpl).format(**d)


class RVT:
    def __init__(self, rvh):
        self._rvh = rvh
        self.var_cmds = {
            "pr": {},
            "rainfall": {},
            "prsn": {},
            "tasmin": {},
            "tasmax": {},
            "tas": {},
            "evspsbl": {},
            "water_volume_transport_in_river_channel": {},
        }
        self.latitude = None
        self.longitude = None
        self.elevation = None
        self.nc_index = 0
        self.grid_weights = None
        self.number_grid_cells = 0

    def hydrate(self, nc_data):
        for var, cmd in nc_data["var_cmds"].items():
            if isinstance(self.var_cmds[var], dict):
                self.var_cmds[var] = replace(cmd, **self.var_cmds[var])
            else:
                self.var_cmds[var] = cmd
        if nc_data["latitude"]:
            self.latitude = nc_data["latitude"]
        else:
            self.latitude = [self._rvh.hrus[0].latitude]
        if nc_data["longitude"]:
            self.longitude = nc_data["longitude"]
        else:
            self.longitude = [self._rvh.hrus[0].longitude]
        if nc_data["elevation"]:
            self.elevation = nc_data["elevation"]
        else:
            self.elevation = [self._rvh.hrus[0].elevation]
        self.number_grid_cells = nc_data["number_grid_cells"]

    def update(self, key, value):
        if key in self.var_cmds:
            self.var_cmds[key].update(value)
            return True
        if key == "nc_index":
            self.nc_index = value
        return False

    def to_rv(self):
        tpl = Path("/home/christian/rpy/ravenpy/models/global/global.rvt").read_text()

        d = {
            "gauge": "",
            "forcing_list": "",
            "observed_data": "",
        }

        use_gauge = any(type(cmd) is DataCommand for cmd in self.var_cmds.values())
        if use_gauge:
            data = []
            for var, cmd in self.var_cmds.items():
                if cmd and var != "water_volume_transport_in_river_channel":
                    data.append(cmd)
            d["gauge"] = GaugeCommand(
                latitude=self.latitude[self.nc_index],
                longitude=self.longitude[self.nc_index],
                elevation=self.elevation[self.nc_index],
                data=data,
            )
        else:
            gw = self.grid_weights or GridWeightsCommand(
                number_hrus=1,
                number_grid_cells=self.number_grid_cells,
                data=[(1, self.nc_index, 1.0)],
            )
            cmds = []
            for var, cmd in self.var_cmds.items():
                if cmd and var != "water_volume_transport_in_river_channel":
                    cmd.grid_weights = gw
                    cmds.append(cmd)
            d["forcing_list"] = "\n".join(map(str, cmds))

        if self.var_cmds.get("water_volume_transport_in_river_channel"):
            d["observed_data"] = self.var_cmds[
                "water_volume_transport_in_river_channel"
            ]

        return tpl.format(**d)


class Config:
    def __init__(self, hrus, subbasins):
        self.rvh = RVH(hrus, subbasins)
        self.rvt = RVT(self.rvh)

    def hydrate(self, rv_type, nc_data):
        self.rvt.hydrate(nc_data)

    def update(self, key, value):
        for rv_type in ("rvh", "rvt"):
            if getattr(self, rv_type).update(key, value):
                break
        return False

    # def to_rv(self):
    #     return self.rvt.to_rv()
