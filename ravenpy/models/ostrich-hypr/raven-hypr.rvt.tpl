#########################################################################
:FileType          rvt ASCII Raven 3.0.4
:WrittenBy         Mahkameh Taheri, Juliane Mai & James Craig
:CreationDate      Feb 2021
#
# Emulation of the HYPR model for Salmon River near Prince George
#------------------------------------------------------------------------

:Gauge meteorological forcings
   :Latitude    {latitude}
   :Longitude   {longitude}
   :Elevation   {elevation}

   :RainCorrection         1.0       #
   :SnowCorrection         par_x21   # para_x21

   {monthly_ave_evaporation}
   {monthly_ave_temperature}

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
