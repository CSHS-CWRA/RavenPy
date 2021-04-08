#########################################################################
:FileType          rvc ASCII Raven 3.0.4
:WrittenBy         Juliane Mai & James Craig
:CreationDate      Feb 2021
#
# Emulation of Blended model simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x29" and "par_x30" wouldn't be detectable)
#    para_half_x29 = para_x29 * 1000. / 2. = par_x29 / 2. [m] = half_x29 [mm]
#    para_half_x30 = para_x30 * 1000. / 2. = par_x30 / 2. [m] = half_x30 [mm]

# initialize to 1/2 full
#:UniformInitialConditions SOIL[0] half_x29 # x(29)*1000/2 [mm]
#:UniformInitialConditions SOIL[1] half_x30 # x(30)*1000/2 [mm]

:HRUStateVariableTable (formerly :IntialConditionsTable)
   :Attributes SOIL[0] SOIL[1]
   :Units mm mm
   1 half_x29 half_x30
:EndHRUStateVariableTable
