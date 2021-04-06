#########################################################################
:FileType          rvt ASCII Raven 3.0.4
:WrittenBy         Robert Chlumsky, James Craig & Juliane Mai
:CreationDate      Feb 2021
#
# Emulation of Canadian Shield simulation of Salmon River near Prince George
#------------------------------------------------------------------------

:Gauge meteorological forcings
   :Latitude    {latitude}
   :Longitude   {longitude}
   :Elevation   {elevation}
   :RainCorrection par_x32  # para_x32
   :SnowCorrection par_x33  # para_x33
   # {raincorrection_cmd}
   # {snowcorrection_cmd}

   {pr}
   {rainfall}
   {prsn}
   {tasmin}
   {tasmax}
   {tas}
   {evspsbl}
:EndGauge

# observed streamflow
{water_volume_transport_in_river_channel}
