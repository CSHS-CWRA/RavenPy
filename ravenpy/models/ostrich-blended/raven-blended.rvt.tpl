#########################################################################
:FileType          rvt ASCII Raven 2.8.2
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Sep 2018
#
# Generic input file for global models with area averaged series
#---------------------------------------------------------------

:Gauge meteorological forcings
   :Latitude    {latitude}
   :Longitude   {longitude}
   :Elevation   {elevation}
   :RainCorrection par_x33
   :SnowCorrection par_x34
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
