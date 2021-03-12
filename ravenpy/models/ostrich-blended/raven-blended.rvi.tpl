#########################################################################
:FileType          rvi ASCII Raven 3.0.4
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Feb 2021
#
# Emulation of Blended model simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#
:Calendar              {calendar}
:RunName               {run_name}-{run_index}
:StartDate             {start_date}
:EndDate               {end_date}
:TimeStep              {time_step}
:Method                ORDERED_SERIES

:PotentialMeltMethod     POTMELT_HMETS
:RainSnowFraction        {rain_snow_fraction}  # RAINSNOW_HBV
:Evaporation             {evaporation}         # PET_OUDIN
:CatchmentRoute          ROUTE_DUMP
:Routing                 ROUTE_NONE
:SoilModel               SOIL_MULTILAYER 3

:Alias DELAYED_RUNOFF CONVOLUTION[1]

:HydrologicProcesses
  :Precipitation   RAVEN_DEFAULT                         ATMOS_PRECIP   MULTIPLE
  :ProcessGroup #infiltration group
                :Infiltration    INF_HMETS               PONDED_WATER   MULTIPLE
                :Infiltration    INF_VIC_ARNO            PONDED_WATER   MULTIPLE
                :Infiltration    INF_HBV                 PONDED_WATER   MULTIPLE
  :EndProcessGroup CALCULATE_WTS par_r01 par_r02
                  :Overflow      OVERFLOW_RAVEN          SOIL[0]        DELAYED_RUNOFF
  :ProcessGroup #quickflow group
                :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[0]        SURFACE_WATER   # interflow, really
                :Baseflow        BASE_VIC                SOIL[0]        SURFACE_WATER
                :Baseflow        BASE_TOPMODEL           SOIL[0]        SURFACE_WATER
  :EndProcessGroup CALCULATE_WTS par_r03 par_r04
  :Percolation                   PERC_LINEAR             SOIL[0]        SOIL[1]         # recharge
    :Overflow                    OVERFLOW_RAVEN          SOIL[1]        DELAYED_RUNOFF
  :Percolation                   PERC_LINEAR             SOIL[1]        SOIL[2]         # loss to deep groundwater (simplifies to HMETS when PERC_COEFF DEEP_GW=0)
  :ProcessGroup #evaporation group
                :SoilEvaporation SOILEVAP_ALL            SOIL[0]        ATMOSPHERE      # AET
                :SoilEvaporation SOILEVAP_TOPMODEL       SOIL[0]        ATMOSPHERE      # AET
  :EndProcessGroup CALCULATE_WTS par_r05
  :Convolve                      CONVOL_GAMMA            CONVOLUTION[0] SURFACE_WATER   # 'surface runoff'
  :Convolve                      CONVOL_GAMMA_2          DELAYED_RUNOFF SURFACE_WATER   # 'delayed runoff'
  :ProcessGroup #quickflow group
                :Baseflow        BASE_LINEAR_ANALYTIC    SOIL[1]        SURFACE_WATER
                :Baseflow        BASE_POWER_LAW          SOIL[1]        SURFACE_WATER
  :EndProcessGroup CALCULATE_WTS par_r06
  :ProcessGroup #snow balance group
                :SnowBalance     SNOBAL_HMETS            MULTIPLE       MULTIPLE
                :SnowBalance     SNOBAL_SIMPLE_MELT      SNOW           PONDED_WATER
                :SnowBalance     SNOBAL_HBV              MULTIPLE       MULTIPLE
                #:SnowBalance    SNOBAL_GAWSER           MULTIPLE       MULTIPLE
  :EndProcessGroup CALCULATE_WTS par_r07 par_r08
:EndHydrologicProcesses

#:CreateRVPTemplate

#---------------------------------------------------------
# Output Options
#
:EvaluationMetrics NASH_SUTCLIFFE RMSE KLING_GUPTA
:SuppressOutput
:SilentMode
:DontWriteWatershedStorage
#
