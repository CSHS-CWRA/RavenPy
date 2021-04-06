#########################################################################
:FileType          rvc ASCII Raven 3.0.4
:WrittenBy         Robert Chlumsky, James Craig & Juliane Mai
:CreationDate      Feb 2021
#
# Emulation of Canadian Shield simulation of Salmon River near Prince George
#------------------------------------------------------------------------
#

# tied parameters:
# (it is important for OSTRICH to find every parameter place holder somewhere in this file)
# (without this "par_x01", "par_x02", and "par_x03" wouldn't be detectable)
#    para_half_x01 = par_half_x01 [mm] = para_x01 [m] * 1000. / 2. = par_x01 * 1000. / 2.
#    para_half_x02 = par_half_x02 [mm] = para_x02 [m] * 1000. / 2. = par_x02 * 1000. / 2.
#    para_half_x03 = par_half_x03 [mm] = para_x03 [m] * 1000. / 2. = par_x03 * 1000. / 2.


:HRUStateVariableTable # (formerly :IntialConditionsTable)
   :Attributes 	       SOIL[0]       SOIL[1]       SOIL[2]
        :Units 	            mm            mm            mm
             1 	  par_half_x01  par_half_x02  par_half_x03
             2 	  par_half_x01          0.00  par_half_x03
#            1   para_half_x01 para_half_x02 para_half_x03
#            2   para_half_x01          0.00 para_half_x03
:EndHRUStateVariableTable