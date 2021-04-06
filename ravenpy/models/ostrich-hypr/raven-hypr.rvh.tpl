#########################################################################
:FileType          rvh ASCII Raven 3.0.4
:WrittenBy         Mahkameh Taheri, Juliane Mai & James Craig
:CreationDate      Feb 2021
#
# Emulation of the HYPR model for Salmon River near Prince George
#------------------------------------------------------------------------
#
#
:SubBasins
        :Attributes     NAME    DOWNSTREAM_ID   PROFILE   REACH_LENGTH    GAUGED
        :Units          none    none            none      km              none
        1,             {name},   -1,             NONE,     _AUTO,          1
:EndSubBasins

:SubBasinProperties
#                                 x_15,
#                               MAXBAS,
  :Parameters,               TIME_CONC,
  :Units     ,  none,                d,
                   1,          par_x15,
#                             para_x15,
:EndSubBasinProperties

{hrus_cmd}
