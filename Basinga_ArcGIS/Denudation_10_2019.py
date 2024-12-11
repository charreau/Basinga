# -*- coding: utf-8 -*-
# ============================================================================================
# Script which calculates the average denudation rate for individual drainage basin
# ============================================================================================

# PYTHON LIBRARIES
# Import arcpy module:
import arcpy
# Basic mathematical functions
from math import sqrt
from Data_base import *

from scipy import interpolate
import os
# Check out any necessary licenses:
arcpy.CheckOutExtension("spatial")
# Set the overwrite option to on:
arcpy.env.overwriteOutput = True

#np.set_printoptions(threshold=np.inf)  # To print a full array
arcpy.AddMessage(np.version.version)

# Attenuation length (g.cm-2)(Braucher et al., EPSL 2011):
lght_att_n = 160.0
lght_att_sm = 1500.0
lght_att_fm = 4320.0
# lght_att_muon = 4814.0  # Mean attenuation length from Braucher, 2013, table 2

# Uncertainties on the 10Be nucleogenic production rates (after Braucher et al. 2013):
U_sm = 0.2
U_fm = 0.2

# Creation of a temporary folder which will contain some intermediate files. They will be deleted
cur_dir = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(cur_dir + "/temp_folder"):
    os.makedirs(cur_dir + "/temp_folder")
    arcpy.AddMessage("temp_folder created")
Output_DEM_file = cur_dir + "/temp_folder/temp_file"

# Script arguments
Input_shape_file_for_erosion_calculations = arcpy.GetParameterAsText(0)
density = float(arcpy.GetParameterAsText(1))  # g/cm3
ConcField = arcpy.GetParameterAsText(2)
C_Uncert = arcpy.GetParameterAsText(3)
If_paleogeomagnetic_correction_previously_applied = arcpy.GetParameter(4)

#=============================================================================================================
# AUXILIARY FUNCTIONS
#=============================================================================================================

def nuclideSelected(nuclide,scaling):
#Table 7Martin et al. (2017) : cREP
    if scaling == 1:#Lal/St scaling + ERA40 atmosphere + Muscheler:
        if nuclide == 10:
            U_spal= 0.048
        elif nuclide == 3: 
            U_spal= 0.0103
        elif nuclide == 21:
            U_spal=0.048
    elif scaling == 0:##LSD scaling + ERA40 atmosphere + Muscheler:
        if nuclide == 10:
            U_spal= 0.063
        elif nuclide == 3:
            U_spal= 0.048
        elif nuclide == 21:
            U_spal=0.063
            
    return U_spal

def optionsSelected(code):
    if code == 1:
        list_fields = ['Pneu_GTI','Pfm_GTI','Psm_GTI']        
    elif code ==2:
        list_fields = ['Pneu_GI','Pfm_GI','Psm_GI']
    elif code ==3:
        list_fields = ['Pneu_TI','Pfm_TI','Psm_TI']
    elif code ==4:
        list_fields = ['Pneu_I','Pfm_I','Psm_I']
    elif code ==5:
        list_fields = ['Pneu_GT','Pfm_GT','Psm_GT']
    elif code == 6:
        list_fields = ['Pneu_G','Pfm_G','Psm_G']
    elif code == 7:
        list_fields = ['Pneu_T','Pfm_T','Psm_T']
    else:
        list_fields = ['Pneu_','Pfm_','Psm_']

    return list_fields

def optionsSelected_VDM(code):
    if code == 1:
        list_fields = ['Pneu_GTP','Pfm_GTI','Psm_GTI']        
    elif code ==2:
        list_fields = ['Pneu_GP','Pfm_GI','Psm_GI']
    elif code ==3:
        list_fields = ['Pneu_TP','Pfm_TI','Psm_TI']
    elif code ==4:
        list_fields = ['Pneu_IP','Pfm_I','Psm_I']
    elif code ==5:
        list_fields = ['Pneu_GTP','Pfm_GT','Psm_GT']
    elif code == 6:
        list_fields = ['Pneu_GP','Pfm_G','Psm_G']
    elif code == 7:
        list_fields = ['Pneu_TP','Pfm_T','Psm_T']
    else:
        list_fields = ['Pneu_P','Pfm_','Psm_']

    return list_fields

#=============================================================================================================
# MAIN
#=============================================================================================================

if If_paleogeomagnetic_correction_previously_applied:
    suffix = "P"
else:
    suffix = ""

# Process: add new fields for erosion rate calculations
arcpy.AddField_management(Input_shape_file_for_erosion_calculations, "denuda_" + suffix, "FLOAT", "", "", 50)
arcpy.AddField_management(Input_shape_file_for_erosion_calculations, "Uncdenu_"  + suffix, "FLOAT", "", "", 50)

rows = arcpy.UpdateCursor(Input_shape_file_for_erosion_calculations)

for row in rows:

    #retrieve the scaling, options selected and nuclide selected during the production rate caculation
    codenuclide=row.getValue("Nuclide")
    codecorrection=row.getValue("Options")
    codescaling=row.getValue("Scaling")

    if If_paleogeomagnetic_correction_previously_applied:
        liste=optionsSelected_VDM(codecorrection)
    else:
        liste=optionsSelected(codecorrection)
    
    # Get the production rates accodringly to the options selected
    ProdR_spal = row.getValue(liste[0])
    ProdR_fm = row.getValue(liste[1])
    ProdR_sm = row.getValue(liste[2])

    #Calculates the uncertainties accordingly to the scaling and nuclide used and based on table 7 of Martin et al. (2007)
    UncProdR_spal = ProdR_spal * nuclideSelected(codenuclide,codescaling)
    UncProdR_sm = ProdR_sm * U_sm
    UncProdR_fm = ProdR_fm * U_fm

    # Get the concentration value and the associated uncertainty
    Conc = float(row.getValue(ConcField))
    U_Conc = float(row.getValue(C_Uncert))

    # Denudation rate calculation
    P_123 = (ProdR_spal*lght_att_n + ProdR_sm*lght_att_sm + ProdR_fm*lght_att_fm) / density
    erosion = P_123 / Conc  # cm/yr

    # Uncertainty calculation
    dc_c2 = U_Conc / (Conc ** 2)
    A = (dc_c2 * P_123) ** 2
    B = (UncProdR_spal * lght_att_n / (Conc * density)) ** 2
    C = (UncProdR_sm * lght_att_sm / (Conc * density)) ** 2
    D = (UncProdR_fm * lght_att_fm / (Conc * density)) ** 2
    Uncert_ero = sqrt(A + B + C + D)

    row.setValue("denuda_"+ suffix, erosion)
    row.setValue("Uncdenu_" + suffix, Uncert_ero)

    rows.updateRow(row)

# If the script is executed while the attribute table is open, fields are not automatically updated.
# There is no working arcpy function to reload the fields, this is a trick to make the refresh automatic.
# A new field is created and immediately removed to force this refresh.
arcpy.AddField_management(Input_shape_file_for_erosion_calculations, "ReloadRows", "TEXT")
arcpy.DeleteField_management(Input_shape_file_for_erosion_calculations, "ReloadRows")
