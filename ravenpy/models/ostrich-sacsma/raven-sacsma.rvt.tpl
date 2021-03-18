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
   :RainCorrection par_x20
   :SnowCorrection par_x21


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
