# -*- coding: utf-8 -*-
# =============================================================================================
# Script which calculates the average cosmogenic production rate for individual drainage basin.
# Based on Stone (2000) model.
# =============================================================================================

# PYTHON LIBRARIES
# Import arcpy module:
import arcpy
# Tools for scientific calculation
from scipy import interpolate
# Basic mathematical functions
from math import *
# Stuff to manage files/folders
import os
import shutil
# Several constants arrays. Numpy is imported with the Data_base file
cur_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(cur_dir)
from Data_base import *
# LSD
from lsd import *

# Check out any necessary licenses:
arcpy.CheckOutExtension("spatial")
# Set the overwrite option to on:
arcpy.env.overwriteOutput = True

# Creation of a temporary folder which will contain some intermediate files. They will be deleted
if not os.path.exists(cur_dir + "/temp_folder"):
    os.makedirs(cur_dir + "/temp_folder")

# =========================================================
#production rate between the different particles at SLHL:
fsp = 0.9886 #Braucher et al.
fmu_slowmuon=0.0027# Braucher et al.
fmu_fastmuon=0.0087# Braucher et al.

# Atmospheric attenuation (g.cm-2)(Braucher et al., EPSL 2011):
at_attenuation_fastmuon = 510
at_attenuation_slowmuon = 260

# Script arguments
Input_shape_file_of_the_drainage_basins = arcpy.GetParameterAsText(0)
Input_DEM_raster = arcpy.GetParameterAsText(1)
Input_nuclide= arcpy.GetParameterAsText(2)
Input_scaling_scheme = arcpy.GetParameterAsText(3)
If_LSD_polynom_check = arcpy.GetParameterAsText(4)
If_geol_mask_checked = arcpy.GetParameter(5)
Input_geol_shape_file = arcpy.GetParameterAsText(6)
For_quartz_zone_extraction = arcpy.GetParameterAsText(7)
Output_geolquartz_polygons2 = arcpy.GetParameterAsText(8)
If_topo_mask_checked = arcpy.GetParameter(9)
Input_topo_Mask_raster = arcpy.GetParameterAsText(10)
If_glims_mask_checked = arcpy.GetParameter(11)
Input_glims_shape_file = arcpy.GetParameterAsText(12)
Output_glims_polygons2 = arcpy.GetParameterAsText(13)
If_VDM_checked = arcpy.GetParameter(14)
#VDM = arcpy.GetParameterAsText(14)
NamFldConc = arcpy.GetParameterAsText(15)
Output_raster_DEM_of_the_drainage_basins = cur_dir + "/temp_folder/temp_file"

VDM = "MUSCHELER"#The code was originally written to propose three possible VDM database but this is now desactivate.
#This can be enabled again by removing this line and uncomment the line "VDM = arcpy.Getparameter...". Be carreful, if yes the properties of the tool must be changed!"

# Definition of the SLHLP : Sea Level High Latitude Production rate (at/g/a):
##10Be:
if Input_nuclide == "10Be":
    if Input_scaling_scheme == "Stone":
        ##Lal/Stone + ERA40 atmosphere + Muscheler:
        SLHLP_Be=4.18 #Table 7 Martin et al. (2017) : cREP

    elif Input_scaling_scheme == "LSD":
        ##LSD + ERA40 atmosphere + Muscheler:
        SLHLP_Be=4.14 #Table 7 Martin et al. (2017) : cREP
        nuclide=10

elif Input_nuclide == "3He":
##Please note that the code was originally written for 10Be only, for simplicty we keep the same for the parameters with "_Be" but it does consider different nuclide
    if Input_scaling_scheme == "Stone":
        ##Lal/Stone + ERA40 atmosphere + Muscheler:
        SLHLP_Be=122 #Table 7Martin et al. (2017) : cREP 

    elif Input_scaling_scheme == "LSD":
        ##LSD + ERA40 atmosphere + Muscheler:
        SLHLP_Be=125 #Table 7 Martin et al. (2017) : cREP
        nuclide=3
        
elif Input_nuclide == "21Ne":
##Please note that the code was originally written for 10Be only, for simplicty we keep the same for the parameters with "_Be" but it does consider different nuclide
    if Input_scaling_scheme == "Stone":
        ##Lal/Stone + ERA40 atmosphere + Muscheler:
        SLHLP_Be=4.18/4.12 #10Be/21Ne=4.12 +/- 0.17 see Kober et al., EPSL 2011 and Table 7 Martin et al. (2017) : cREP

    elif Input_scaling_scheme == "LSD":
        ##LSD + ERA40 atmosphere + Muscheler:
        SLHLP_Be=4.14/4.12 #10Be/21Ne=4.12 +/- 0.17 see Kober et al., EPSL 2011 and Table 7 Martin et al. (2017) : cREP
        nuclide=10
        
SLHLP_Be_neutron = fsp*SLHLP_Be
SLHLP_Be_fastmuon=fmu_fastmuon*SLHLP_Be
SLHLP_Be_slowmuon=fmu_slowmuon*SLHLP_Be

#===========================================================================================================================================================================================
#Auxiliary functions needed for the calculations
#==========================================================================================================================================================================================
def Rc_correctedStoneFactors(VDM, LAT_average, P_average, Time):
    """
    Returns the time-corrected Stone factors
    Inputs:
        VDM         : virtual dipole moment
        LAT_average : average latitude - degrees
        P_average   : average pressure - hPa
        Time        : sample age - years
        VDM is a matrix
        LAT_average, P_average and Time are scalars
    Outputs:
        FactStoneCorrected   : corrected total Stone factor
        FactStoneCorrected_n : corrected neutron Stone factor
        FactStoneCorrected_s : corrected slow muon Stone factor
        FactStoneCorrected_f : corrected fast muon Stone factor
    """
    Dates, Paleo = VDM[0], VDM[1]
    sea_level_P = 1013.25
    LAT_average_rad = LAT_average * math.pi / 180

    StoneFact_n = []
    
    for p in Paleo:
        Rc = (p * 1e22 * 4 * 1e-7 * 3 * 1e8) / (16 * 1e9 * (6.3712 * 1e6) ** 2) * (np.cos(LAT_average_rad)) ** 4
        stone_n = Rc_StoneFactors(Rc, P_average, sea_level_P)
        StoneFact_n.append(stone_n)

    # Time average of Stone' factors by integration
    StoneFactAve_n = AveragIntegr(Dates, StoneFact_n)

    StoneFactAve_n[0] = StoneFactAve_n[1]

    #Interpolation near the value
    NbPoints = 500
    #Mid = np.floor(NbPoints / 2)
    Mid=250
    VectT = np.linspace(Time * 0.5, Time * 1.5, NbPoints)

    VectSt_n = interpolate.splrep(Dates, StoneFactAve_n)

    VectS_n = interpolate.splev(VectT, VectSt_n)

    VecTF = [VectT[Mid - 1], VectT[Mid], VectT[Mid + 1]]
    VecTEnd = []

    for TF in VecTF:
        VecTEnd.append(abs(TF - Time))

    minimum = min(VecTEnd)
    Indice = VecTEnd.index(minimum)
    Indice = Mid + Indice - 2

    FactStoneCorrected_n = VectS_n[Indice]
    arcpy.AddMessage(FactStoneCorrected_n)  

    return  FactStoneCorrected_n

def Rc_StoneFactors(Rc, P, SLP):
    """
    Interpolates the Stone factors corresponding to Rc
    Inputs:
        Rc  : cutoff rigidity
        P   : pressure - hPa
        SLP : sea Level Pressure - hPa
        Rc, P and SLP are scalars
    Outputs:
        StoneFact_total    : total Stone factor
        StoneFact_neutron  : neutron Stone factor
        StoneFact_slowmuon : slow muon Stone factor
        StoneFact_fastmuon : fast muon Stone factor
    Note:
        Modified from CREp code - StonefactorL.m
    """
    # Latitudes in degrees and Stone coefficients
    i_lat_d, a, b, c, d, e = St[:, 0], St[:, 1], St[:, 2], St[:, 3], St[:, 4], St[:, 5]

    # Latitudes for interpolation, in radians
    i_lat_r = i_lat_d * np.pi / 180

    # Cutoff rigidities for interpolation
    i_Rc = 14.3 * np.cos(i_lat_r) ** 4

    # Scaling factors for interpolation
    i_sf = a + b * np.exp(-P / 150) + c * P + d * P ** 2 + e * P ** 3

    # Extending to 0 rigidity
    i_Rc = np.append(i_Rc, 0)
    i_sf = np.append(i_sf, i_sf[6])

    # Extending to 23GV by fitting a log-log line to 0-20 latitude values, i.e. Rc > 10GV
    fits = np.polyfit(np.log(i_Rc[:3]), np.log(i_sf[:3]), 1)
    add_sf = np.exp(np.log(i_sf[0]) + fits[0] * (np.log(np.arange(25, 15, -1)) - np.log(i_Rc[0])))
    i_sf = np.insert(i_sf, 0, add_sf)
    add_rc = np.arange(25, 15, -1)
    i_Rc = np.insert(i_Rc, 0, add_rc)

    # Use of splrep requires to have i_Rc in ascending order. Then, i_sf needs to be reversed
    i_Rc = np.sort(i_Rc)
    i_sf = i_sf[::-1]

    # Interpolation
    spline_representation = interpolate.splrep(i_Rc, i_sf)
    StoneFact_neutron = interpolate.splev(Rc, spline_representation)

    return  StoneFact_neutron

def AveragIntegr(Dates, StoneFact):
    """
    Time integration using trapezoidal rule
    """
    #arcpy.AddMessage(StoneFact)
    YIntegr = []
    YAveVect = []
    for i, date in enumerate(Dates):
        DeltaX = date - Dates[i - 1]
        Trapeze = 0.5 * (StoneFact[i] + StoneFact[i - 1]) 

        y_int = YIntegr[i - 1] + DeltaX * Trapeze if i != 0 else 0
        y_ave = y_int / (date - Dates[0]) if i != 0 else 0

        if i == 0:
            y_int = StoneFact[0] 
        
        YIntegr.append(y_int)
        YAveVect.append(y_ave)

    YAveVect[0]=StoneFact[0] 
    Corr_fact=YAveVect/StoneFact[0]   
    return Corr_fact

def processVDM(input_file, VDM_in, input_correction_type):
    """
    Computes the paleomagnetic correction
    Inputs:
        input_file            : .shp file to work with
        VDM_in                : Virtual Dipole Moment used
        input_correction_type : correction applied, it is a string containing either "G", "I", "T" or any conbination.
    """
    arcpy.AddMessage("Processing paleomagnetic correction")
    
    if VDM_in == "GLOPIS":
        VDM_in = GLOPIS
        VDM_in_m=GLOPISm
    elif VDM_in == "MUSCHELER":
        VDM_in = MUSCHELER
        VDM_in_m=MUSCHELERm
    elif VDM_in == "VDM2000":
        VDM_in = VDM2000
        VDM_in_m=VDM2000m

    sfn_non_corr='SFn_' + input_correction_type

    suffix = input_correction_type + "P"  # Example: "GIP"
    

    if Input_scaling_scheme == "LSD":
        addFields(input_file, ['Time', 'Pneu_'+ suffix, 'Psm_'+ suffix, 'Pfm_'+ suffix,  'SFn_' + suffix, 'SFf_' + suffix, 'SFs_' + suffix])
    else:
        addFields(input_file, ['Time', 'Pneu_'+ suffix,'SFn_' + suffix])
        
    Teq = "[" + NamFldConc + "]/[Pneu_" + input_correction_type + "]"

    arcpy.CalculateField_management(input_file, "Time", Teq, "VB", "")

    rows = arcpy.UpdateCursor(input_file)

    for row in rows:
        arcpy.AddMessage("Processing row " + str(row.getValue("FID")))
 
        latitude = row.getValue("LAT_ave")
        longitude = row.getValue("LONG_ave")
        Zaverage=row.getValue("Alt_ave")
        time = row.getValue("Time")
        Sfn_non_corr=row.getValue(sfn_non_corr)

        # Computes corrected Stone factors

        if Input_scaling_scheme == "LSD":
            SFn_corr,SFmu_corr=timecorrect_LSD(latitude,longitude,Zaverage,time,VDM_in_m)
            
            Pneu=SFn_corr*SLHLP_Be_neutron
            Psm=SFmu_corr*SLHLP_Be_slowmuon
            Pfm=SFmu_corr*SLHLP_Be_fastmuon
            Pneu=Pneu[0]
            Pfm=Pfm[0]
            Psm=Psm[0]
            
            row.setValue("Pneu_"+ suffix, Pneu)
            row.setValue("Psm_"+ suffix, Psm)
            row.setValue("Pfm_"+ suffix, Pfm)
            row.setValue("SFn_" + suffix, SFn_corr[0])
            row.setValue("SFs_" + suffix, SFmu_corr[0])
            row.setValue("SFf_" + suffix, SFmu_corr[0])        
            
        else:
            pressure = row.getValue("Pres")
            pressure =1013.25
            Corr_factor = Rc_correctedStoneFactors(VDM_in, latitude, pressure, time)#gives a correction factor not a stone factor...
            
            SFn_corr=Sfn_non_corr*Corr_factor

            Pneu=SFn_corr*SLHLP_Be_neutron
            
            row.setValue("Pneu_"+ suffix, Pneu)
            row.setValue("SFn_" + suffix, SFn_corr)
            
        rows.updateRow(row)

def StoneCoefficients(latitude):
    """
    Takes a latitude in parameter and returns the associated Stone coefficients.
    Inputs:
        latitude : latitude where the Stone coefficient are wanted - degrees
    Outputs:
        a, b, c, d, e : Stone coefficients
    Note: 
        M = (dec * St[0, 6] + (10 - dec) * St[1, 6]) / 10 is not used because Muons scaling factors are described following Braucher et al.(2011)
        The absolute value of latitude is taken, so a south hemisphere area is treated as if it was a north area.
    """
    latitude = abs(latitude)
    if latitude < 10:
        dec = 10 - latitude
        a = (dec * St[0, 1] + (10 - dec) * St[1, 1]) / 10
        b = (dec * St[0, 2] + (10 - dec) * St[1, 2]) / 10
        c = (dec * St[0, 3] + (10 - dec) * St[1, 3]) / 10
        d = (dec * St[0, 4] + (10 - dec) * St[1, 4]) / 10
        e = (dec * St[0, 5] + (10 - dec) * St[1, 5]) / 10
    elif latitude < 20:
        dec = 20 - latitude
        a = (dec * St[1, 1] + (10 - dec) * St[2, 1]) / 10
        b = (dec * St[1, 2] + (10 - dec) * St[2, 2]) / 10
        c = (dec * St[1, 3] + (10 - dec) * St[2, 3]) / 10
        d = (dec * St[1, 4] + (10 - dec) * St[2, 4]) / 10
        e = (dec * St[1, 5] + (10 - dec) * St[2, 5]) / 10
    elif latitude < 30:
        dec = 30 - latitude
        a = (dec * St[2, 1] + (10 - dec) * St[3, 1]) / 10
        b = (dec * St[2, 2] + (10 - dec) * St[3, 2]) / 10
        c = (dec * St[2, 3] + (10 - dec) * St[3, 3]) / 10
        d = (dec * St[2, 4] + (10 - dec) * St[3, 4]) / 10
        e = (dec * St[2, 5] + (10 - dec) * St[3, 5]) / 10
    elif latitude < 40:
        dec = 40 - latitude
        a = (dec * St[3, 1] + (10 - dec) * St[4, 1]) / 10
        b = (dec * St[3, 2] + (10 - dec) * St[4, 2]) / 10
        c = (dec * St[3, 3] + (10 - dec) * St[4, 3]) / 10
        d = (dec * St[3, 4] + (10 - dec) * St[4, 4]) / 10
        e = (dec * St[3, 5] + (10 - dec) * St[4, 5]) / 10
    elif latitude < 50:
        dec = 50 - latitude
        a = (dec * St[4, 1] + (10 - dec) * St[5, 1]) / 10
        b = (dec * St[4, 2] + (10 - dec) * St[5, 2]) / 10
        c = (dec * St[4, 3] + (10 - dec) * St[5, 3]) / 10
        d = (dec * St[4, 4] + (10 - dec) * St[5, 4]) / 10
        e = (dec * St[4, 5] + (10 - dec) * St[5, 5]) / 10
    elif latitude < 60:
        dec = 60 - latitude
        a = (dec * St[5, 1] + (10 - dec) * St[6, 1]) / 10
        b = (dec * St[5, 2] + (10 - dec) * St[6, 2]) / 10
        c = (dec * St[5, 3] + (10 - dec) * St[6, 3]) / 10
        d = (dec * St[5, 4] + (10 - dec) * St[6, 4]) / 10
        e = (dec * St[5, 5] + (10 - dec) * St[6, 5]) / 10
    else:  # LAT >= 60
        a = St[6, 1]
        b = St[6, 2]
        c = St[6, 3]
        d = St[6, 4]
        e = St[6, 5]

    return a, b, c, d, e


def ERA40_Data(latitude, longitude):
    """
    Returns the corresponding sea level pressure and temperature of a site from the ERA40 database, and the adiabatic lapse rate
    Inputs:
        latitude  : latitude of the site - degrees
        Longitude : longitude of the site - degrees
    Outputs:
        sea_level_P : sea level pressure for given latitude and longitude from the ERA40 database - hPa
        sea_level_T : sea level temperature for given latitude and longitude from the ERA40 database - K
        dtdz        : Lifton Adiabatic Lapse Rate Fit to COSPAR CIRA-86 <10 km altitude - K/m
    """
    sea_level_P = float(interpolate.interp2d(ERA40_lon, ERA40_lat, ERA40_meanP).__call__(longitude, latitude))
    sea_level_T = float(interpolate.interp2d(ERA40_lon, ERA40_lat, ERA40_meanT).__call__(longitude, latitude))
    dtdz = LIFTON_lr[0] + LIFTON_lr[1] * latitude + LIFTON_lr[2] * latitude ** 2 + LIFTON_lr[3] * latitude ** 3 \
           + LIFTON_lr[4] * latitude ** 4 + LIFTON_lr[5] * latitude ** 5 + LIFTON_lr[6] * latitude ** 6
    dtdz = -dtdz

    return sea_level_P, sea_level_T, dtdz


def DEMtoStoneCoefficients(DEM_array, bot_lat, cellsize):
    """
    Computes the Stone coefficients for each cell.
    Inputs:
        DEM_array : the numpy array corresponding to the raster
        bot_lat   : the south latitude of the raster - degrees
        cellsize  : the cellsize of the raster
    Ouputs:
        a, b, c, d, e : 1D-matrices of all the Stone coefficients.
    Note:
        NoData-cells are not returned
    """
    a = DEM_array.copy()
    b = DEM_array.copy()
    c = DEM_array.copy()
    d = DEM_array.copy()
    e = DEM_array.copy()
    for row in range(0, DEM_array.shape[0]):
        a[row], b[row], c[row], d[row], e[row] = StoneCoefficients(bot_lat + row * cellsize)

    return a[DEM_array > 0], b[DEM_array > 0], c[DEM_array > 0], d[DEM_array > 0], e[DEM_array > 0]

def correction_latitude (latitude,cellsize):

    """
    Computes the surface weight for latitudinal correction

    Input: latitude

    Output: weight

    """
    Rt=6370000 #Earth radius in m
    
    non_zero = latitude != 0
    lat_min=np.min(abs(latitude[non_zero]))
    lat_rad=latitude[non_zero]*pi/180
    r_min=Rt*np.cos(np.radians(lat_min))
    r=Rt*np.cos(lat_rad)
    Circ_min=2*pi*r_min
    Circ=2*pi*r
    Circ_lon=2*pi*Rt
    d_1deg_lon=Circ_lon/360
    d_1deg_lat=Circ/360
    d_1deg_lat_min=Circ_min/360

    cellsize_lon=cellsize*d_1deg_lon
    cellsize_min_lat=cellsize*d_1deg_lat_min
    surf_ref=cellsize_lon*cellsize_min_lat
    
    cellsize_all=cellsize*d_1deg_lat

    Surf_all=cellsize_all*cellsize_lon

    rate=np.abs(Surf_all/surf_ref)

    return rate

def timecorrect_LSD(lat_ave,long_ave,Z_ave,AGE,VDM):
    """
    Correct basin average the scaling factor from LSD for paleomagnetic variation

    Input:
        conc: basin average concentration
        VDM: paleomagnetic data_base
        10Be: basin average cosmogenic production rate
        lat_ave: basin average latitude
        long_ave: basin average longitude
        Z_ave: basin average elevation
    Outputs:
        ScalingFactorIntegre: scaling factor integrated during time
    """
    #tmp_SF = LSDtimeInteg(10, tmp_time_vec, tmp_SF_vec)
    time_vec, scalingFactors_vecSP,scalingFactors_vecMU = LSDv9(lat_ave, long_ave, Z_ave, 0, AGE, -1, nuclide, VDM)
    SFn_integre = LSDtimeInteg(10, time_vec, scalingFactors_vecSP)
    SFmu_integre = LSDtimeInteg(10, time_vec, scalingFactors_vecMU)
    
    return SFn_integre,SFmu_integre

def getRasterProperties(DEM):
    """
    Extracts some properties from a DEM file
    Input:
        DEM : raster file
    Outputs:
        ncols      : number of columns in the input raster
        nrows      : number of rows in the input raster
        bot_lat    : south latitude - degrees
        left_long  : west longitude - degrees
        top_lat    : north latitude - degrees
        right_long : east longitude - degrees
        cellsize   : x-size of pixels
    """        
    ncols      = float( arcpy.GetRasterProperties_management(DEM, "COLUMNCOUNT").getOutput(0) )
    nrows      = float( arcpy.GetRasterProperties_management(DEM, "ROWCOUNT"   ).getOutput(0) )
    bot_lat    = float( arcpy.GetRasterProperties_management(DEM, "BOTTOM"     ).getOutput(0) )
    left_long  = float( arcpy.GetRasterProperties_management(DEM, "LEFT"       ).getOutput(0) )
    top_lat    = float( arcpy.GetRasterProperties_management(DEM, "TOP"        ).getOutput(0) )
    right_long = float( arcpy.GetRasterProperties_management(DEM, "RIGHT"      ).getOutput(0) )
    cellsize   = float( arcpy.GetRasterProperties_management(DEM, "CELLSIZEX"  ).getOutput(0) )

    return ncols, nrows, bot_lat, left_long, top_lat, right_long, cellsize


def addFields(shape_file, fields_list):
    """
    Adds the fields contained in fields_list to the input shape_file, as a FLOAT of length 50*
    Inputs:
        shape_file  : .shp file
        fields_list : list of the fields to add to the shapefile
    """
    for field in fields_list:
        arcpy.AddField_management(shape_file, field, "FLOAT", "", "", 50)


def removeFields(shape_file, fields_list):
    """
    Removes all the fields contained in fields_list from the input shape_file
    Inputs:
        shape_file  : .shp file
        fields_list : list of the fields to delete from the shapefile
    """
    for field in fields_list:
        arcpy.DeleteField_management(shape_file, field)

def pixel_count(input_basins_shp, input_DEM,output_DEM):
    """
    Computes the number of pixels for GLIMS correction.
    Inputs:
        input_basins_shp : .shp file with basins
        input_DEM        : raster file
        output_DEM       : temporary rasters which will be deleted
    Output:
        input_basins_shp : .shp file with basins
    """
    # Local variables:
    clipmask_shp = 'this_name_doesnt_matter'
    FID = "FID"

    # Process: add new fields for production rate results and associated parameters to the attribute table of input_basins_shp
    
    addFields(input_basins_shp, ['pixels'])

    # Read the number of lines/features in the file input_basins_shp exported above
    rows = arcpy.UpdateCursor(input_basins_shp)

    # Loop over each line to export raster of DEM
    for row in rows:
        a = row.FID
        
        # Set the output raster file
        output_DEM_ite = output_DEM
        # Set the FID number
        where_clause = '"' + FID + '"' + '= %d' % a
        # Select features according to their FID number
        arcpy.Select_analysis(input_basins_shp, clipmask_shp, where_clause)
        # Process: Extraction by DEM
        arcpy.gp.ExtractByMask_sa(input_DEM, clipmask_shp, output_DEM_ite)

        DEM_array = arcpy.RasterToNumPyArray(output_DEM_ite, nodata_to_value=-9999)
        DEM_array = DEM_array.astype(float)

        pixel_with_data = DEM_array >= 0  # Boolean array that stores where there is Data pixel

        pixels = np.sum(DEM_array > 0)
        row.setValue("pixels", pixels)

        rows.updateRow(row)

    return input_basins_shp

def code_type(option):
    """
    Function to return a code for the option that is stored in the attribute table
    """
    arcpy.AddMessage(option)
    if option == 'GTI':
        code_option=1
    elif option == 'GI':
        code_option=2
    elif option == 'TI':
        code_option=3
    elif option == 'I':
        code_option=4
    elif option == 'GT':
        code_option=5
    elif option == 'G':
        code_option=6
    elif option == 'T':
        code_option=7
    elif option == 'None' :
        code_option=8

    code_option=float(code_option)
    return code_option

def code_nuclide(Imput_nuclide):
    """
    Function to return a code for nuclide that is stored in the attribute table
    """
    if Imput_nuclide == '10Be':
        code_nucl=10
    elif Imput_nuclide == '3He':
        code_nucl=3
    elif Imput_nuclide == '21Ne':
        code_nucl=21
    
    code_nucl=float(code_nucl)
    return code_nucl

def code_scaling(Input_scaling_scheme):
    """
    Function to return a code for sclaing that is stored in the attribute table
    """   
    if Input_scaling_scheme == 'ls':
        code_scal=1
    else:
        code_scal=0
        
    code_scal=float(code_scal)    
    return code_scal

#===========================================================================================================================================================================================
#Function to calculate the scaling factors and production rates accordingly either to the Lal/St or to the LSD model
#=========================================================================================================================================================================================== 

def processNumpy(DEM, MASK=None):
    """
    Computes the sum of several parameters
    Inputs:
        DEM            : raster file
        MASK(optional) : topographic mask
    Outputs:
            SUM_Be   : total production rates - g/at/yr
            SUM_Ben  : neutrons production rates - g/at/yr
            SUM_SFn  : neutrons Stone factors
            SUM_Besm : slow muons production rates - g/at/yr
            SUM_SFsm : slow muons Stone factors
            SUM_Befm : fast muons production rates - g/at/yr
            SUM_SFfm : fast muons Stone factors
            SUM_P    : pressures - hPa
            SUM_LAT  : latitudes - degrees
            SUM_LONG : longitudes - degrees
            SUM_DEM  : altitudes of the DEM - m
            pixels   : number of pixels with data
            nrows    : number of rows
            SUM_MASK : topographic factors
    Note:
        - All the outputs are sums over all pixels.
        - If a pixel from the DEM is negative, it will be seen as a NoData value and will not be counted. This may lead to unexpected errors.
        - The data pixels from the mask has to fit exactly with the data pixels from the DEM.
    """
    DEM_array = arcpy.RasterToNumPyArray(DEM, nodata_to_value=-9999)
    DEM_array = DEM_array.astype(float)

    pixel_with_data = DEM_array >= 0  # Boolean array that stores where there is Data pixel

    taille0=np.shape(DEM_array)

    if MASK is not None:  # Test if there is a mask
        MASK_array = arcpy.RasterToNumPyArray(MASK, nodata_to_value=-9999)
        MASK_array = MASK_array.astype(float)

        if not np.array_equal(pixel_with_data, MASK_array >= 0):  # Test if both NoData-Data pixels are the same between DEM and MASK
            arcpy.AddError("<!> DEM and MASK don't have the exact same pixels with data <!>\n"
                            + "DEM:" + str(DEM_array[pixel_with_data].size) + " vs MASK:" + str(MASK_array[MASK_array >= 0].size)
                            + "The results may be wrong.")

    else:  # There is no MASK, MASK_array is the same shape than the DEM_array, filled with ones
        MASK_array = np.ones(DEM_array.shape)

    ncols, nrows, bot_lat, left_long, top_lat, right_long, cellsize = getRasterProperties(DEM)

    if taille0[0] != nrows:#correction made for bug in the function RasterToNumpy which sometime may generate additional rows when compared to the true original extracted DEM
        nrows=taille0[0]

    # Sea level pressure, sea level temperature and adiabatic lapse rate at the center of the ASCII file from the ERA40 atmosphere model
    center_LAT = bot_lat + 0.5 * nrows * cellsize
    center_LONG = left_long + 0.5 * ncols * cellsize
    sea_level_P, sea_level_T, dtdz = ERA40_Data(center_LAT, center_LONG)

    a, b, c, d, e = DEMtoStoneCoefficients(DEM_array, bot_lat, cellsize)
    
    LAT = np.linspace(bot_lat, top_lat, nrows).reshape(nrows, 1) * pixel_with_data  # Latitudes array with zeros where DEM is negative

    rate=correction_latitude(LAT,cellsize)#weighting for latitudinal variation of pixel surface
    #arcpy.AddMessage(np.mean(rate))

    # PRESSURE TABLE as a function of elevation (Stone, 2000):
    # Nodata_values are excluded
    P = sea_level_P * np.exp(-0.03417 / dtdz * (np.log(sea_level_T) - np.log(sea_level_T - dtdz * DEM_array[pixel_with_data])))

    SF_Be_neutron = (a + b * np.exp(-P / 150) + c * P + d * P ** 2 + e * P ** 3) * MASK_array[pixel_with_data]
    SF_Be_neutron=SF_Be_neutron*rate
    ProdDEM_Be_neutron = SLHLP_Be_neutron * SF_Be_neutron

    SF_Be_fastmuon = np.exp((1013.25 - P) / at_attenuation_fastmuon) * MASK_array[pixel_with_data]
    SF_Be_fastmuon = SF_Be_fastmuon *rate
    ProdDEM_Be_fastmuon = SLHLP_Be_fastmuon * SF_Be_fastmuon

    SF_Be_slowmuon = np.exp((1013.25 - P) / at_attenuation_slowmuon) * MASK_array[pixel_with_data]
    SF_Be_slowmuon = SF_Be_slowmuon *rate
    ProdDEM_Be_slowmuon = SLHLP_Be_slowmuon * SF_Be_slowmuon

    pixels = np.sum(DEM_array > 0)
    
    SUM_LAT = np.sum(LAT)
    LONG = np.linspace(left_long, right_long, ncols).reshape(1, ncols) * pixel_with_data  # Latitudes array with zeros where DEM is negative
    SUM_LONG = np.sum(LONG)
    SUM_DEM = np.sum(DEM_array[pixel_with_data])
    SUM_Pn = np.sum(ProdDEM_Be_neutron)
    SUM_Psm = np.sum(ProdDEM_Be_slowmuon)
    SUM_Pfm = np.sum(ProdDEM_Be_fastmuon)
    SUM_SFn = np.sum(SF_Be_neutron)
    SUM_SFsm = np.sum(SF_Be_slowmuon)
    SUM_SFfm = np.sum(SF_Be_fastmuon)
    SUM_MASK = np.sum(MASK_array[pixel_with_data])
    SUM_Pre=np.sum(P)

    if pixels == 0:
        arcpy.AddError("<!> No pixel with data found in the basin <!>")

    return pixels, SUM_Pre, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK

def processNumpy_LSD_poly(DEM, MASK=None):
    """
    LSD scaling factor using a polynomial interpolation form Stone factors
    Computes the sum of several parameters
    Inputs:
        conc           : name of the concentration field
        DEM            : raster file
        MASK(optional) : topographic mask
    Outputs:
        SUM_Be   : total production rates following the LSD scaling model- g/at/yr
        SUM_LAT  : latitudes - degrees
        SUM_LONG : longitudes - degrees
        SUM_DEM  : altitudes of the DEM - m
        pixels   : number of pixels with data
        nrows    : number of rows
        SUM_MASK : topographic factors
    Note:
        It seems that the LSD code doesn't take muons into account. As a consequence, SUM_Be might be equivalent to SUM_Ben if using Stone.
    """
    import random
    DEM_array = arcpy.RasterToNumPyArray(DEM, nodata_to_value=-9999)
    DEM_array = DEM_array.astype(float)
    
    pixel_with_data = DEM_array >= 0  # Boolean array that stores where there is Data pixel

    DEM_bon=DEM_array[pixel_with_data]
  
    taille0=np.shape(DEM_array)

    if MASK is not None:  # Test if there is a mask
        MASK_array = arcpy.RasterToNumPyArray(MASK, nodata_to_value=-9999)
        MASK_array = MASK_array.astype(float)

        if not np.array_equal(pixel_with_data, MASK_array >= 0):  # Test if both NoData-Data pixels are the same between DEM and MASK
            arcpy.AddError("<!> DEM and MASK don't have the exact same pixels with data <!>\n"
                            + "DEM:" + str(DEM_array[pixel_with_data].size) + " vs MASK:" + str(MASK_array[MASK_array >= 0].size)
                            + "The results may be wrong.")

    else:  # There is no MASK, MASK_array is the same shape than the DEM_array, filled with ones
        MASK_array = np.ones(DEM_array.shape)

    ncols, nrows, bot_lat, left_long, top_lat, right_long, cellsize = getRasterProperties(DEM)

    if taille0[0] != nrows:#correction made for bug in the function RasterToNumpy which sometime may generate additional rows when compared to the true original extracted DEM
        nrows=taille0[0]

    LAT = np.linspace(bot_lat, top_lat, nrows).reshape(nrows, 1) * pixel_with_data  # Latitudes array with zeros where DEM is negative
    LONG = np.linspace(left_long, right_long, ncols).reshape(1, ncols) * pixel_with_data  # Longitudes array with zeros where DEM is negative

    LAT_bon=LAT[pixel_with_data]
    LONG_bon=LONG[pixel_with_data]

    rate=correction_latitude(LAT,cellsize)#weighting for latitudinal variation of pixel surface
    
    ##Step 1: calculate Stone factors on all pixels:
    # Sea level pressure, sea level temperature and adiabatic lapse rate at the center of the ASCII file from the ERA40 atmosphere model
    center_LAT = bot_lat + 0.5 * nrows * cellsize
    center_LONG = left_long + 0.5 * ncols * cellsize
    sea_level_P, sea_level_T, dtdz = ERA40_Data(center_LAT, center_LONG)

    a, b, c, d, e = DEMtoStoneCoefficients(DEM_array, bot_lat, cellsize)
    
    # PRESSURE TABLE as a function of elevation (Stone, 2000):
    # Nodata_values are excluded
    P = sea_level_P * np.exp(-0.03417 / dtdz * (np.log(sea_level_T) - np.log(sea_level_T - dtdz * DEM_array[pixel_with_data])))
    
    SF_Be_neutron = (a + b * np.exp(-P / 150) + c * P + d * P ** 2 + e * P ** 3)# * MASK_array[pixel_with_data]

    SF_Be_fastmuon = np.exp((1013.25 - P) / at_attenuation_fastmuon) #* MASK_array[pixel_with_data]

    SF_Be_slowmuon = np.exp((1013.25 - P) / at_attenuation_slowmuon) #* MASK_array[pixel_with_data]

    Total_S_st = (SF_Be_neutron*fsp) + (SF_Be_fastmuon*fmu_fastmuon) + (SF_Be_slowmuon*fmu_slowmuon)
    
    ##Step 2:random selection of DEM data to calculate the LSD factors:
    
    length=len(Total_S_st)
    
    if length >20000:
        n_echantillon=1000
    else :
        n_echantillon=np.ceil(length/10)  

    interval=(max(Total_S_st)-min(Total_S_st))/10
    mini=min(Total_S_st)
    
    St_sp_extract_end=[]
    St_sm_extract_end=[]
    St_fm_extract_end=[]
    DEM_extract_end=[]
    Lat_extract_end=[]
    Lon_extract_end=[]
    
    for o in range(1,11):
        indy=(Total_S_st>mini+((o-1)*interval)) & (Total_S_st<mini+(o*interval))  
        
        if len(Total_S_st[indy])>n_echantillon/10:
            St_sp_extract=random.sample(SF_Be_neutron[indy],n_echantillon/10)
            St_sm_extract=random.sample(SF_Be_slowmuon[indy],n_echantillon/10)
            St_fm_extract=random.sample(SF_Be_fastmuon[indy],n_echantillon/10)
            DEM_extract=random.sample(DEM_bon[indy],n_echantillon/10)
            Lat_extract=random.sample(LAT_bon[indy],n_echantillon/10)
            Lon_extract=random.sample(LONG_bon[indy],n_echantillon/10)

        else:
            St_sp_extract=SF_Be_neutron[indy]
            St_sm_extract=SF_Be_slowmuon[indy]
            St_fm_extract=SF_Be_fastmuon[indy]
            DEM_extract=DEM_bon[indy]
            Lat_extract=LAT_bon[indy]
            Lon_extract=LONG_bon[indy]

        St_sp_extract_end= np.concatenate((St_sp_extract_end,St_sp_extract))
        St_sm_extract_end= np.concatenate((St_sm_extract_end,St_sm_extract))
        St_fm_extract_end= np.concatenate((St_fm_extract_end,St_fm_extract))
        DEM_extract_end= np.concatenate((DEM_extract_end,DEM_extract))
        Lat_extract_end= np.concatenate((Lat_extract_end,Lat_extract))
        Lon_extract_end= np.concatenate((Lon_extract_end,Lon_extract))

    ##Step 3: calculate the LSD factors on the sampled cells of the DEM:
    
    lsd_sp_extract_end = []
    lsd_mu_extract_end = []
    
    arcpy.AddMessage("Calculate LSD on sampled cells")
    arcpy.AddMessage("Number of sampled cells:")
    arcpy.AddMessage(len(St_sp_extract_end))

    for i, dem in enumerate(DEM_extract_end):  # Loop over the DEM array, because LSD code is not fully vectorized and needs single values as inputs
        if dem <= 0: continue # Skip DEM negative values
            
        tmp_time_vec, tmp_SF_vecSP, tmp_SF_vecMU= LSDv9(Lat_extract_end[i], Lon_extract_end[i], dem, 0, 0, -1, nuclide, MUSCHELERm)
     
        lsd_sp_extract_end.append(tmp_SF_vecSP[0])#*MASK_array.flatten()[i])
        lsd_mu_extract_end.append(tmp_SF_vecMU[0])#*MASK_array.flatten()[i])

    ##Step 4: polynomial fit between Stone and LSD factors on the sampled pixels:

    arcpy.AddMessage("Calculate the polyfit")
    Preg = np.polyfit(St_sp_extract_end, lsd_sp_extract_end, 4)# determine the polynom to calculte SP LSD production from stone production
    Preg2 = np.polyfit(St_sm_extract_end, lsd_mu_extract_end, 4)# determine the polynom to calculte slow muons LSD production from Braucher production
    Preg3 = np.polyfit(St_fm_extract_end, lsd_mu_extract_end, 4)# determine the polynom to calculte fast muons LSD production from Braucher production

    ##Step 5:extrapollation of  LSD factors to all cells based on the polifit relation:
    #Total_Prod_LSD = []
    #for i in range(len(Total_Prod_st)):
    #    Total_Prod_LSD.append(Preg[0] * Total_Prod_st[i] ** 4 + Preg[1] * Total_Prod_st[i] ** 3 + Preg[2] * Total_Prod_st[i] ** 2 + Preg[3] * Total_Prod_st[i] + Preg[4])
    arcpy.AddMessage("Extrapolate to all cells")
    Total_LSD_sp=Preg[0] * SF_Be_neutron ** 4 + Preg[1] * SF_Be_neutron ** 3 + Preg[2] * SF_Be_neutron ** 2 + Preg[3] * SF_Be_neutron + Preg[4]
    Total_LSD_sp=Total_LSD_sp*rate* MASK_array[pixel_with_data]
    Total_LSD_sm=Preg2[0] * SF_Be_slowmuon ** 4 + Preg2[1] * SF_Be_slowmuon ** 3 + Preg2[2] * SF_Be_slowmuon ** 2 + Preg2[3] * SF_Be_slowmuon + Preg2[4]
    Total_LSD_sm=Total_LSD_sm*rate* MASK_array[pixel_with_data]
    Total_LSD_fm=Preg3[0] * SF_Be_fastmuon ** 4 + Preg3[1] * SF_Be_fastmuon ** 3 + Preg3[2] * SF_Be_fastmuon ** 2 + Preg3[3] * SF_Be_fastmuon + Preg3[4]
    Total_LSD_fm=Total_LSD_fm*rate* MASK_array[pixel_with_data]

    ProdDEM_Be_neutron= Total_LSD_sp*SLHLP_Be_neutron
    ProdDEM_Be_neutron=ProdDEM_Be_neutron*rate
    ProdDEM_Be_slowmuon= Total_LSD_sm*SLHLP_Be_slowmuon
    ProdDEM_Be_slowmuon=ProdDEM_Be_slowmuon*rate
    ProdDEM_Be_fastmuon=Total_LSD_fm*SLHLP_Be_fastmuon
    ProdDEM_Be_fastmuon=ProdDEM_Be_fastmuon*rate
   
    pixels = np.sum(DEM_array > 0)
    SUM_LAT = np.sum(LAT)
    SUM_LONG = np.sum(LONG)
    SUM_DEM = np.sum(DEM_array[pixel_with_data])
    SUM_MASK = np.sum(MASK_array[pixel_with_data])
    SUM_Ben = np.sum(ProdDEM_Be_neutron)
    SUM_Besm = np.sum(ProdDEM_Be_slowmuon)
    SUM_Befm = np.sum(ProdDEM_Be_fastmuon)
    SUM_SFn=np.sum(Total_LSD_sp)
    SUM_SFsm=np.sum(Total_LSD_sm)
    SUM_SFfm=np.sum(Total_LSD_fm)

    return pixels,SUM_LAT, SUM_LONG, SUM_DEM,SUM_Ben, SUM_Besm, SUM_Befm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK

def processNumpy_LSD(DEM, MASK=None):
    """
    Raw LSD sclaing factors
    Computes the sum of several parameters
    Inputs:
        conc           : name of the concentration field
        DEM            : raster file
        MASK(optional) : topographic mask
    Outputs:
        SUM_Be   : total production rates following the LSD scaling model- g/at/yr
        SUM_LAT  : latitudes - degrees
        SUM_LONG : longitudes - degrees
        SUM_DEM  : altitudes of the DEM - m
        pixels   : number of pixels with data
        nrows    : number of rows
        SUM_MASK : topographic factors
    Note:
        It seems that the LSD code doesn't take muons into account. As a consequence, SUM_Be might be equivalent to SUM_Ben if using Stone.
    """
    DEM_array = arcpy.RasterToNumPyArray(DEM, nodata_to_value=-9999)
    DEM_array = DEM_array.astype(float)

    pixel_with_data = DEM_array >= 0  # Boolean array that stores where there is Data pixel

    taille0=np.shape(DEM_array)

    if MASK is not None:  # Test if there is a mask
        MASK_array = arcpy.RasterToNumPyArray(MASK, nodata_to_value=-9999)
        MASK_array = MASK_array.astype(float)

        if not np.array_equal(pixel_with_data, MASK_array >= 0):  # Test if both NoData-Data pixels are the same between DEM and MASK
            arcpy.AddError("<!> DEM and MASK don't have the exact same pixels with data <!>\n"
                            + "DEM:" + str(DEM_array[pixel_with_data].size) + " vs MASK:" + str(MASK_array[MASK_array >= 0].size)
                            + "The results may be wrong.")

    else:  # There is no MASK, MASK_array is the same shape than the DEM_array, filled with ones
        MASK_array = np.ones(DEM_array.shape)

    ncols, nrows, bot_lat, left_long, top_lat, right_long, cellsize = getRasterProperties(DEM)

    if taille0[0] != nrows:#Correction made for bug in the function RasterToNumpy which sometime may generate additional rows when compared to the true original extracted DEM
        nrows=taille0[0]

    LAT = np.linspace(bot_lat, top_lat, nrows).reshape(nrows, 1) * pixel_with_data  # Latitudes array with zeros where DEM is negative
    LONG = np.linspace(left_long, right_long, ncols).reshape(1, ncols) * pixel_with_data  # Longitudes array with zeros where DEM is negative

    SUM_Pn = 0
    SUM_Psm = 0
    SUM_Pfm = 0
    SUM_SFn = 0
    SUM_SFsm = 0
    SUM_SFfm = 0

    rate=correction_latitude(LAT,cellsize)#weighting for latitudinal variation of pixel surface

    for i, dem in enumerate(np.nditer(DEM_array[pixel_with_data])):  # Loop over the DEM array, because LSD code is not fully vectorized and needs single values as inputs
  
        tmp_time_vec, tmp_SF_vecSP, tmp_SF_vecMU = LSDv9(LAT.flatten()[i], LONG.flatten()[i], dem, 0, 0, -1, nuclide, MUSCHELERm)
        tmp_SF_vecSP=tmp_SF_vecSP*rate[i]
        tmp_SF_vecMU=tmp_SF_vecMU*rate[i]
        
        SUM_SFn +=tmp_SF_vecSP[0]*MASK_array.flatten()[i] 
        SUM_SFsm +=tmp_SF_vecMU[0]*MASK_array.flatten()[i]
        SUM_SFfm +=tmp_SF_vecMU[0]*MASK_array.flatten()[i]
        SUM_Pn +=tmp_SF_vecSP[0] *SLHLP_Be_neutron
        SUM_Psm +=tmp_SF_vecMU[0] *SLHLP_Be_fastmuon
        SUM_Pfm +=tmp_SF_vecMU[0] *SLHLP_Be_slowmuon

    pixels = np.sum(DEM_array > 0)
    SUM_LAT = np.sum(LAT)
    SUM_LONG = np.sum(LONG)
    SUM_DEM = np.sum(DEM_array[pixel_with_data])
    SUM_MASK = np.sum(MASK_array[pixel_with_data])

    return pixels, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK

#============================================================================================================================================================================================================
#Main functions used to extract the D.E.M, the elevations, latitudes and calculate the production rates and scaling factors accordingly to the optionS selected (TS maks, geological mask, glacial coer, VDM)
#============================================================================================================================================================================================================

def geol_topo(input_basins_shp, input_DEM, input_geol_shp, input_query, input_MASK, output_DEM, output_geol,input_correction_type,Input_nuclide,Input_scaling_scheme):
    """
    Computes the production rate with geological and topographic masks.
    Inputs:
        input_basins_shp : .shp file with basins
        input_DEM        : raster file
        input_geol_shp   : .shp with geological mask
        input_query      : geological query
        input_MASK       : topographic mask
        output_DEM       : temporary rasters which will be deleted
        output_geol      : .shp that will be created with the geological mask
    Output:
        output_shp : new .shp file
    """
    # Local variables:
    clipgeol_shp = 'this_name_doesnt_matter'

    # Process: If necessary, select in the appropriate field the geological soils of quartz nature using an SQL expression
    arcpy.Select_analysis(input_geol_shp, clipgeol_shp, input_query)

    # Process: Create the new quartz polygon shape file inside each drainage basin
    arcpy.Clip_analysis(input_basins_shp, clipgeol_shp, output_geol, "")
    output_shp = output_geol + '.shp'

    # ===========================================================================================================================================
    # Process: Iterate the selection of all individual feature (i.e drainage basins) of the new quartz polygon shape file of the drainage basins.
    # The goal is to extract the DEM and topographic mask of each individual basin in order to further calculate the production rates.
    # ================================================================================================================================

    # Local variables:
    clipmask_shp = 'this_name_doesnt_matter'
    FID = "FID"

    # Process: add new fields for production rate results and associated parameters to the attribute table of the Output_geolquartz_polygons
    
    if Input_scaling_scheme == "Stone":
        addFields(output_shp, ['pres','LAT_ave', 'LONG_ave','Alt_ave','Pneu_'+input_correction_type, 'Psm_'+input_correction_type, 'Pfm_'+input_correction_type,'SFn_'+input_correction_type, 'SFs_'+input_correction_type, 'SFf_'+input_correction_type,'Topo_fc','Nuclide','Scaling','Options'])
    elif Input_scaling_scheme == "LSD":
        addFields(output_shp, ['LAT_ave', 'LONG_ave','Alt_ave','Pneu_'+input_correction_type, 'Psm_'+input_correction_type, 'Pfm_'+input_correction_type,'SFn_'+input_correction_type, 'SFs_'+input_correction_type, 'SFf_'+input_correction_type,'Topo_fc','Nuclide','Scaling','Options'])
  
    # Read the number of lines/features in the file 'Output_geolquartz_polygons' exported above
    rows = arcpy.UpdateCursor(output_shp)

    codenucl = code_nuclide(Input_nuclide)  
    codescal = code_scaling(Input_scaling_scheme)
    optioncode = code_type(input_correction_type)

    # Loop over each line to export raster of DEM and topographic mask
    for row in rows:
        a = row.FID
        arcpy.AddMessage("Processing row " + str(a))

        # Set the output raster file
        output_DEM_ite = output_DEM
        output_MASK_ite = output_DEM + '_ts'  # ts: for topographic shielding
        # Set the FID number
        where_clause = '"' + FID + '"' + '= %d' % a
        # Select features according to their FID number
        arcpy.Select_analysis(output_shp, clipmask_shp, where_clause)
        # Process: Extraction by DEM and topographic masks
        arcpy.gp.ExtractByMask_sa(input_DEM, clipmask_shp, output_DEM_ite)
        arcpy.gp.ExtractByMask_sa(input_MASK, clipmask_shp, output_MASK_ite)

        if Input_scaling_scheme == "Stone":
            # Process calculations on pixels
            pp, pressure, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm ,SUM_MASK \
                = processNumpy(output_DEM_ite,output_MASK_ite)  # No mask, we don't need the last value returned

            # Calculation of average production rates in the drainage basins
            row.setValue("LAT_ave", SUM_LAT / pp)
            row.setValue("LONG_ave", SUM_LONG / pp)
            row.setValue("Alt_ave", SUM_DEM / pp)
            row.setValue("Pneu_"+input_correction_type, SUM_Pn / pp)
            row.setValue("Psm_"+input_correction_type, SUM_Psm / pp)
            row.setValue("Pfm_"+input_correction_type, SUM_Pfm / pp)
            row.setValue("SFn_"+input_correction_type, SUM_SFn / pp)
            row.setValue("SFs_"+input_correction_type, SUM_SFsm / pp)
            row.setValue("SFf_"+input_correction_type, SUM_SFfm / pp)
            row.setValue("Topo_fc", SUM_MASK / pp)
            row.setValue("Pres", pressure / pp)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        elif Input_scaling_scheme == "LSD":
            # Process calculations on pixels
            if If_LSD_polynom_check:
                arcpy.AddMessage("LSD: Polynomial simplification")
                pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm,SUM_MASK  = processNumpy_LSD_poly(output_DEM_ite,output_MASK_ite)

            else:
                pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn1, SUM_Psm1, SUM_Pfm1, SUM_SFn1, SUM_SFsm1, SUM_SFfm1,SUM_MASK  = processNumpy_LSD(output_DEM_ite,output_MASK_ite)
                SUM_Pn=SUM_Pn1[0]
                SUM_Psm=SUM_Psm1[0]
                SUM_Pfm=SUM_Pfm1[0]
                SUM_SFn=SUM_SFn1[0]
                SUM_SFsm=SUM_SFsm1[0]
                SUM_SFfm=SUM_SFfm1[0]

            row.setValue("LAT_ave", SUM_LAT / pp)
            row.setValue("LONG_ave", SUM_LONG / pp)
            row.setValue("Alt_ave", SUM_DEM / pp)
            row.setValue("Pneu_GT", SUM_Pn/ pp)
            row.setValue("Psm_GT", SUM_Psm/ pp)
            row.setValue("Pfm_GT", SUM_Pfm/ pp)
            row.setValue("SFn_GT", SUM_SFn / pp)
            row.setValue("SFs_GT", SUM_SFsm / pp)
            row.setValue("SFf_GT", SUM_SFfm/ pp)
            row.setValue("Topo_fc", SUM_MASK / pp)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        rows.updateRow(row)

    return output_shp


def geol_notopo(input_basins_shp, input_DEM, input_geol_shp, input_query, output_DEM, output_geol, input_correction_type,Input_nuclide,Input_scaling_scheme):
    """
    Computes the production rate with a geological mask.
    Inputs:
        input_basins_shp : .shp file with basins
        input_DEM        : raster file
        input_geol_shp   : .shp with geological mask
        input_query      : geological query
        output_DEM       : temporary rasters which will be deleted
        output_geol      : .shp that will be created with the geological mask
    Output:
        output_shp : new .shp file
    """
    # Local variables:
    clipgeol_shp = 'this_name_doesnt_matter'

    # Process: If necessary, select in the appropriate field the geological soils of quartz nature using an SQL expression
    arcpy.Select_analysis(input_geol_shp, clipgeol_shp, input_query)

    # Process: Create the new quartz polygon shape file inside each drainage basin
    arcpy.Clip_analysis(input_basins_shp, clipgeol_shp, output_geol, "")
    output_shp = output_geol + '.shp'

    # =============================================================================================================================
    # Process: Iterate the selection of all individual feature (i.e drainage basins) of the quartz polygons of the drainage basins.
    # The goal is to extract the DEM of each individual basin in order to further calculate the production rates.
    # ===========================================================================================================

    # Local variables:
    clipmask_shp = "this_name_doesnt_matter"
    FID = "FID"

    # Process: add new fields for production rate results and associated parameters to the attribute table of the Output_geolquartz_polygons
    if Input_scaling_scheme == "Stone":
        addFields(output_shp, ['Pres','LAT_ave', 'LONG_ave','Alt_ave','Pres','Pneu_'+input_correction_type, 'Psm_'+input_correction_type, 'Pfm_'+input_correction_type,'SFn_'+input_correction_type, 'SFs_'+input_correction_type, 'SFf_'+input_correction_type,'Nuclide','Scaling','Options'])
    elif Input_scaling_scheme == "LSD":
        addFields(output_shp, ['LAT_ave', 'LONG_ave','Alt_ave','Pres','Pneu_'+input_correction_type, 'Psm_'+input_correction_type, 'Pfm_'+input_correction_type,'SFn_'+input_correction_type, 'SFs_'+input_correction_type, 'SFf_'+input_correction_type,'Nuclide','Scaling','Options'])
    
    # Read the number of lines/features in the file 'Output_geolquartz_polygons' exported above
    rows = arcpy.UpdateCursor(output_shp)

    codenucl = code_nuclide(Input_nuclide)  
    codescal = code_scaling(Input_scaling_scheme)
    optioncode = code_type(input_correction_type)

    # Loop over each line to export raster of DEM
    for row in rows:
        a = row.FID
        arcpy.AddMessage("Processing row " + str(a))

        # Set the output raster file
        output_DEM_ite = output_DEM
        # Set the FID number
        where_clause = '"' + FID + '"' + '= %d' % a
        # Select features according to their FID number
        arcpy.Select_analysis(output_shp, clipmask_shp, where_clause)
        # Process: Extraction by DEM mask
        arcpy.gp.ExtractByMask_sa(input_DEM, clipmask_shp, output_DEM_ite)

        if Input_scaling_scheme == "Stone":
            # Process calculations on pixels
            pp, pressure, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm \
                = processNumpy(output_DEM_ite)[:-1]  # No mask, we don't need the last value returned  # No mask, we don't need the last value returned

            # Calculation of average production rates in the drainage basins
            row.setValue("LAT_ave", SUM_LAT / pp)
            row.setValue("LONG_ave", SUM_LONG / pp)
            row.setValue("Alt_ave", SUM_DEM / pp)
            row.setValue("Pneu_"+input_correction_type, SUM_Pn / pp)
            row.setValue("Psm_"+input_correction_type, SUM_Psm / pp)
            row.setValue("Pfm_"+input_correction_type, SUM_Pfm / pp)
            row.setValue("SFn_"+input_correction_type, SUM_SFn / pp)
            row.setValue("SFs_"+input_correction_type, SUM_SFsm / pp)
            row.setValue("SFf_"+input_correction_type, SUM_SFfm / pp)
            row.setValue("Pres", pressure / pp)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)
            

        elif Input_scaling_scheme == "LSD":
            # Process calculations on pixels
            if If_LSD_polynom_check:
                arcpy.AddMessage("LSD: Polynomial simplification")
                pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm = processNumpy_LSD_poly(output_DEM_ite)[:-1]

            else:
                pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn1, SUM_Psm1, SUM_Pfm1, SUM_SFn1, SUM_SFsm1, SUM_SFfm1 = processNumpy_LSD(output_DEM_ite)[:-1]
                SUM_Pn=SUM_Pn1[0]
                SUM_Psm=SUM_Psm1[0]
                SUM_Pfm=SUM_Pfm1[0]
                SUM_SFn=SUM_SFn1[0]
                SUM_SFsm=SUM_SFsm1[0]
                SUM_SFfm=SUM_SFfm1[0]

            row.setValue("LAT_ave", SUM_LAT / pp)
            row.setValue("LONG_ave", SUM_LONG / pp)
            row.setValue("Alt_ave", SUM_DEM / pp)
            row.setValue("Pneu_"+input_correction_type, SUM_Pn/ pp)
            row.setValue("Psm_"+input_correction_type, SUM_Psm/ pp)
            row.setValue("Pfm_"+input_correction_type, SUM_Pfm/ pp)
            row.setValue("SFn_"+input_correction_type, SUM_SFn / pp)
            row.setValue("SFs_"+input_correction_type, SUM_SFsm / pp)
            row.setValue("SFf_"+input_correction_type, SUM_SFfm/ pp)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        rows.updateRow(row)

    return output_shp


def nogeol_topo(input_basins_shp, input_DEM, input_MASK, output_DEM, input_correction_type,Input_nuclide,Input_scaling_scheme):
    """
    Computes the production rate with a topographic masks.
    Inputs:
        input_basins_shp : .shp file with basins
        input_DEM        : raster file
        input_MASK       : topographic mask
        output_DEM       : temporary rasters which will be deleted
    Output:
        input_basins_shp : .shp file with basins
    """
    # Local variables:
    clipmask_shp = 'this_name_doesnt_matter'
    FID = "FID"

    # Process: add new fields for production rate results and associated parameters to the attribute table of the Input_shape_file_of_the_drainage_basins

    if Input_scaling_scheme == "Stone":
        addFields(input_basins_shp, ['Pres','LAT_ave', 'LONG_ave','Alt_ave','Pneu_'+ input_correction_type, 'Psm_'+ input_correction_type, 'Pfm_'+ input_correction_type,'SFn_'+ input_correction_type, 'SFs_'+ input_correction_type, 'SFf_'+ input_correction_type,'Topo_fc','Nuclide','Scaling','Options'])
    elif Input_scaling_scheme == "LSD":
        addFields(input_basins_shp, ['LAT_ave', 'LONG_ave','Alt_ave','Pneu_'+ input_correction_type, 'Psm_'+ input_correction_type, 'Pfm_'+ input_correction_type,'SFn_'+ input_correction_type, 'SFs_'+ input_correction_type, 'SFf_'+ input_correction_type,'Topo_fc','Nuclide','Scaling','Options'])
       
    # Read the number of lines/features in the file 'Input_shape_file_of_the_drainage_basins' exported above
    rows = arcpy.UpdateCursor(input_basins_shp)
    
    codenucl = code_nuclide(Input_nuclide)  
    codescal = code_scaling(Input_scaling_scheme)
    optioncode = code_type(input_correction_type)

    # Loop over each line to export raster of DEM and topographic mask
    for row in rows:
        a = row.FID
        arcpy.AddMessage("Processing row " + str(a))

        # Set the output raster file
        output_DEM_ite = output_DEM
        output_MASK_ite = output_DEM + '_ts'  # ts: for topographic shielding
        # Set the FID number
        where_clause = '"' + FID + '"' + '= %d' % a
        # Select features according to their FID number
        arcpy.Select_analysis(input_basins_shp, clipmask_shp, where_clause)
        # Process: Extraction by DEM and topographic masks
        arcpy.gp.ExtractByMask_sa(input_DEM, clipmask_shp, output_DEM_ite)
        arcpy.gp.ExtractByMask_sa(input_MASK, clipmask_shp, output_MASK_ite)

        if Input_scaling_scheme == "Stone":
            # Process calculations on pixels
            pp, pressure, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm ,SUM_MASK \
                = processNumpy(output_DEM_ite,output_MASK_ite)  # No mask, we don't need the last value returned

            # Calculation of average production rates in the drainage basins
            row.setValue("LAT_ave", SUM_LAT / pp)
            row.setValue("LONG_ave", SUM_LONG / pp)
            row.setValue("Alt_ave", SUM_DEM / pp)
            row.setValue("Pneu_"+ input_correction_type, SUM_Pn / pp)
            row.setValue("Psm_"+ input_correction_type, SUM_Psm / pp)
            row.setValue("Pfm_"+ input_correction_type, SUM_Pfm / pp)
            row.setValue("SFn_"+ input_correction_type, SUM_SFn / pp)
            row.setValue("SFs_"+ input_correction_type, SUM_SFsm / pp)
            row.setValue("SFf_"+ input_correction_type, SUM_SFfm / pp)
            row.setValue("Topo_fc", SUM_MASK / pp)
            row.setValue("Pres", pressure / pp)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        elif Input_scaling_scheme == "LSD":
            # Process calculations on pixels
            if If_LSD_polynom_check:
                arcpy.AddMessage("LSD: Polynomial simplification")
                pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm,SUM_MASK  = processNumpy_LSD_poly(output_DEM_ite,output_MASK_ite)

            else:
                pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn1, SUM_Psm1, SUM_Pfm1, SUM_SFn1, SUM_SFsm1, SUM_SFfm1,SUM_MASK  = processNumpy_LSD(output_DEM_ite,output_MASK_ite)
                SUM_Pn=SUM_Pn1[0]
                SUM_Psm=SUM_Psm1[0]
                SUM_Pfm=SUM_Pfm1[0]
                SUM_SFn=SUM_SFn1[0]
                SUM_SFsm=SUM_SFsm1[0]
                SUM_SFfm=SUM_SFfm1[0]

            row.setValue("LAT_ave", SUM_LAT / pp)
            row.setValue("LONG_ave", SUM_LONG / pp)
            row.setValue("Alt_ave", SUM_DEM / pp)
            row.setValue("Pneu_"+ input_correction_type, SUM_Pn/ pp)
            row.setValue("Psm_"+ input_correction_type, SUM_Psm/ pp)
            row.setValue("Pfm_"+ input_correction_type, SUM_Pfm/ pp)
            row.setValue("SFn_"+ input_correction_type, SUM_SFn / pp)
            row.setValue("SFs_"+ input_correction_type, SUM_SFsm / pp)
            row.setValue("SFf_"+ input_correction_type, SUM_SFfm/ pp)
            row.setValue("Topo_fc", SUM_MASK / pp)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        rows.updateRow(row)    

    return input_basins_shp

def nogeol_notopo(input_basins_shp, input_DEM, output_DEM,Input_nuclide,Input_scaling_scheme):

    """
    Computes the production rate without any mask.
    Inputs:
        input_basins_shp : .shp file with basins
        input_DEM        : raster file
        output_DEM       : temporary rasters which will be deleted
    Output:
        input_basins_shp : .shp file with basins
    """
    # Local variables:
    clipmask_shp = 'this_name_doesnt_matter'
    FID = "FID"

    # Process: add new fields for production rate results and associated parameters to the attribute table of input_basins_shp
    if Input_scaling_scheme == "Stone":
        addFields(input_basins_shp, ['Pres','LAT_ave', 'LONG_ave','Alt_ave','Pneu_', 'Psm_', 'Pfm_','SFn_', 'SFs_', 'SFf_','Nuclide','Scaling','Options'])
    elif Input_scaling_scheme == "LSD":
        addFields(input_basins_shp, ['LAT_ave', 'LONG_ave','Alt_ave','Pneu_', 'Psm_', 'Pfm_','SFn_', 'SFs_', 'SFf_','Nuclide','Scaling','Options'])
    # Read the number of lines/features in the file input_basins_shp exported above
    rows = arcpy.UpdateCursor(input_basins_shp)

    codenucl = code_nuclide(Input_nuclide)  
    codescal = code_scaling(Input_scaling_scheme)
    optioncode = code_type('None')

    # Loop over each line to export raster of DEM
    for row in rows:
        a = row.FID
        arcpy.AddMessage("Processing row " + str(a))
        
        # Set the output raster file
        output_DEM_ite = output_DEM
        # Set the FID number
        where_clause = '"' + FID + '"' + '= %d' % a
        # Select features according to their FID number
        arcpy.Select_analysis(input_basins_shp, clipmask_shp, where_clause)
        # Process: Extraction by DEM
        arcpy.gp.ExtractByMask_sa(input_DEM, clipmask_shp, output_DEM_ite)

        if Input_scaling_scheme == "Stone":
            # Process calculations on pixels
            pp, pressure, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm \
                = processNumpy(output_DEM_ite)[:-1]  # No mask, we don't need the last value returned

            # Calculation of average production rates in the drainage basins
            row.setValue("LAT_ave", SUM_LAT / pp)
            row.setValue("LONG_ave", SUM_LONG / pp)
            row.setValue("Alt_ave", SUM_DEM / pp)
            row.setValue("Pneu_", SUM_Pn / pp)
            row.setValue("Psm_", SUM_Psm / pp)
            row.setValue("Pfm_", SUM_Pfm / pp)
            row.setValue("SFn_", SUM_SFn / pp)
            row.setValue("SFs_", SUM_SFsm / pp)
            row.setValue("SFf_", SUM_SFfm / pp)
            row.setValue("Pres", pressure / pp)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        elif Input_scaling_scheme == "LSD":
            # Process calculations on pixels
            if If_LSD_polynom_check:
                arcpy.AddMessage("LSD: Polynomial simplification")
                pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm = processNumpy_LSD_poly(output_DEM_ite)[:-1]

            else:
                pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn1, SUM_Psm1, SUM_Pfm1, SUM_SFn1, SUM_SFsm1, SUM_SFfm1 = processNumpy_LSD(output_DEM_ite)[:-1]
                SUM_Pn=SUM_Pn1[0]
                SUM_Psm=SUM_Psm1[0]
                SUM_Pfm=SUM_Pfm1[0]
                SUM_SFn=SUM_SFn1[0]
                SUM_SFsm=SUM_SFsm1[0]
                SUM_SFfm=SUM_SFfm1[0]

            row.setValue("LAT_ave", SUM_LAT / pp)
            row.setValue("LONG_ave", SUM_LONG / pp)
            row.setValue("Alt_ave", SUM_DEM / pp)
            row.setValue("Pneu_", SUM_Pn/ pp)
            row.setValue("Psm_", SUM_Psm/ pp)
            row.setValue("Pfm_", SUM_Pfm/ pp)
            row.setValue("SFn_", SUM_SFn / pp)
            row.setValue("SFs_", SUM_SFsm / pp)
            row.setValue("SFf_", SUM_SFfm/ pp)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        rows.updateRow(row)

    return input_basins_shp


def glims_topo(input_shp, input_DEM, input_MASK, input_glims_shp, output_DEM, output_glims, input_correction_type,Input_nuclide,Input_scaling_scheme):
    """
    Computes the ice cover correction with a topographic mask
    Inputs:
        input_shp             : .shp file with basins
        input_DEM             : raster file
        input_MASK            : topographic mask
        input_glims_shp       : .shp file with the ice cover
        output_DEM            : temporary rasters which will be deleted
        output_glims          : new .shp that will be created with ice correction
        input_correction_type : correction applied, it is a string containing either "G", "I", "T" or any conbination.
    Output:
        output_shp : new .shp file
    """
    
    
    # Process: Create the new glims polygon shape file inside each drainage basin
    arcpy.Erase_analysis(input_shp, input_glims_shp, output_glims, "")
    output_shp = output_glims + '.shp'
    
    # ======================================================================================================================================
    # Process: Iterate the selection of all individual feature (i.e drainage basins) of the ice free quartz polygons of the drainage basins.
    # The goal is to extract the DEM and topographic mask of each individual basin in order to further calculate the production rates.
    # ================================================================================================================================

    # Local variables:
    clipmask_shp = 'this_name_doesnt_matter'
    FID = "FID"

    # Process: add new fields for production rate results and associated parameters to the attribute table of the Output_glims_polygons
  
    if Input_scaling_scheme == "Stone":
        addFields(output_shp, ['Pres','LAT_ave', 'LONG_ave','Alt_ave','Pneu_'+ input_correction_type, 'Psm_'+ input_correction_type, 'Pfm_'+ input_correction_type,'SFn_'+ input_correction_type, 'SFs_'+ input_correction_type, 'SFf_'+ input_correction_type,'Topo_fc','pixfreeI','Nuclide','Scaling','Options'])
    elif Input_scaling_scheme == "LSD":
        addFields(output_shp, ['LAT_ave', 'LONG_ave','Alt_ave','Pneu_'+ input_correction_type, 'Psm_'+ input_correction_type, 'Pfm_'+ input_correction_type,'SFn_'+ input_correction_type, 'SFs_'+ input_correction_type, 'SFf_'+ input_correction_type,'Topo_fc','pixfreeI','Nuclide','Scaling','Options'])

    # Read the number of lines/features in the file 'Output_glims_polygons' exported above
    rows = arcpy.UpdateCursor(output_shp)

    codenucl = code_nuclide(Input_nuclide)  
    codescal = code_scaling(Input_scaling_scheme)
    optioncode = code_type(input_correction_type)

    # Loop over each line to export raster of DEM and topographic mask
    for row in rows:
        a = row.FID
        arcpy.AddMessage("Processing row " + str(a))

        # Set the output raster file
        output_DEM_ite = output_DEM
        output_MASK_ite = output_DEM + '_ts'  # ts: for topographic shielding
        # Set the OBJECTID number
        where_clause = '"' + FID + '"' + '= %d' % a
        # Select features according to their OBJECTID number
        arcpy.Select_analysis(output_shp, clipmask_shp, where_clause)
        # Process: Extraction of the DEM and topographic masks
        arcpy.gp.ExtractByMask_sa(input_DEM, clipmask_shp, output_DEM_ite)
        arcpy.gp.ExtractByMask_sa(input_MASK, clipmask_shp, output_MASK_ite)

        pixels = row.getValue("pixels")

        if Input_scaling_scheme == "Stone":
            # Process calculations on pixels
            ppg, pressure, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK  \
                = processNumpy(output_DEM_ite, output_MASK_ite)  

            # Calculation of average production rates in the drainage basins:
            row.setValue("LAT_ave", SUM_LAT / ppg)
            row.setValue("LONG_ave", SUM_LONG / ppg)
            row.setValue("Alt_ave", SUM_DEM / ppg)
            row.setValue("Pneu_"+ input_correction_type, SUM_Pn / pixels)
            row.setValue("Psm_"+ input_correction_type, SUM_Psm / pixels)
            row.setValue("Pfm_"+ input_correction_type, SUM_Pfm / pixels)
            row.setValue("SFn_" + input_correction_type, SUM_SFn / pixels)
            row.setValue("SFf_" + input_correction_type, SUM_SFfm / pixels)
            row.setValue("SFs_" + input_correction_type, SUM_SFsm / pixels)
            row.setValue("Topo_fc", SUM_MASK / ppg)
            row.setValue("pixfreeI",  ppg)
            row.setValue("Pres", pressure / ppg)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)
            
        elif Input_scaling_scheme == "LSD":
            # Process calculations on pixels
            if If_LSD_polynom_check:
                arcpy.AddMessage("LSD: Polynomial simplification")
                ppg, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK = processNumpy_LSD_poly(output_DEM_ite, output_MASK_ite)
            else:
                ppg, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK  = processNumpy_LSD(output_DEM_ite, output_MASK_ite)
                SUM_Pn=SUM_Pn1[0]
                SUM_Psm=SUM_Psm1[0]
                SUM_Pfm=SUM_Pfm1[0]
                SUM_SFn=SUM_SFn1[0]
                SUM_SFsm=SUM_SFsm1[0]
                SUM_SFfm=SUM_SFfm1[0]
                
            row.setValue("LAT_ave", SUM_LAT / ppg)
            row.setValue("LONG_ave", SUM_LONG / ppg)
            row.setValue("Alt_ave", SUM_DEM / ppg)
            row.setValue("Pneu_"+ input_correction_type, SUM_Pn / pixels)
            row.setValue("Psm_"+ input_correction_type, SUM_Psm / pixels)
            row.setValue("Pfm_"+ input_correction_type, SUM_Pfm / pixels)
            row.setValue("SFn_" + input_correction_type, SUM_SFn / pixels)
            row.setValue("SFf_" + input_correction_type, SUM_SFfm / pixels)
            row.setValue("SFs_" + input_correction_type, SUM_SFsm / pixels)
            row.setValue("Topo_fc", SUM_MASK / ppg)
            row.setValue("pixfreeI",  ppg)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        rows.updateRow(row)

    return output_shp


def glims_notopo(input_shp, input_DEM, input_glims_shp, output_DEM, output_glims, input_correction_type,Input_nuclide,Input_scaling_scheme):
    """
    Computes the ice cover correction without the topographic mask
    Inputs:
        input_shp             : .shp file with basins
        input_DEM             : raster file
        input_glims_shp       : .shp file with the ice cover
        output_DEM            : temporary rasters which will be deleted
        output_glims          : new .shp that will be created with ice correction
        input_correction_type : correction applied, it is a string containing either "G", "I", "T" or any conbination.
    Output:
        output_shp : new .shp file
    """
    # Process: Create the new glims polygon shape file inside each drainage basin
    arcpy.Erase_analysis(input_shp, input_glims_shp, output_glims, "")
    output_shp = output_glims + '.shp'

    # ====================================================================================================================================================
    # Process: Iterate the selection of all individual feature (i.e drainage basins) of the new ice free quartz polygon shape file of the drainage basins.
    # The goal is to extract the DEM of each individual basin in order to further calculate the production rates.
    # ===========================================================================================================

    # Local variables:
    clipmask_shp = 'this_name_doesnt_matter'
    FID = "FID"

    # Process: add new fields for production rate results and associated parameters to the attribute table of the Output_glims_polygons
    
    if Input_scaling_scheme == "Stone":
        addFields(output_shp, ['pres','LAT_ave', 'LONG_ave','Alt_ave','Pneu_'+ input_correction_type, 'Psm_'+ input_correction_type, 'Pfm_'+ input_correction_type,'SFn_'+ input_correction_type, 'SFs_'+ input_correction_type, 'SFf_'+ input_correction_type,'pixfreeI','Nuclide','Scaling','Options'])
    elif Input_scaling_scheme == "LSD":
        addFields(output_shp, ['LAT_ave', 'LONG_ave','Alt_ave','Pneu_'+ input_correction_type, 'Psm_'+ input_correction_type, 'Pfm_'+ input_correction_type,'SFn_'+ input_correction_type, 'SFs_'+ input_correction_type, 'SFf_'+ input_correction_type,'pixfreeI','Nuclide','Scaling','Options'])
 
    # Read the number of lines/features in the file 'Output_glims_polygons' exported above
    rows = arcpy.UpdateCursor(output_shp)

    codenucl = code_nuclide(Input_nuclide)  
    codescal = code_scaling(Input_scaling_scheme)
    optioncode = code_type(input_correction_type)

    # Loop over each line to export raster of DEM
    for row in rows:
        a = row.FID
        arcpy.AddMessage("Processing row " + str(a))

        # Set the output raster file
        output_DEM_ite = output_DEM + '%d' % a
        # Set the OBJECTID number
        where_clause = '"' + FID + '"' + '= %d' % a
        # Select features according to their OBJECTID number
        arcpy.Select_analysis(output_shp, clipmask_shp, where_clause)
        # Process: Extraction by DEM mask
        arcpy.gp.ExtractByMask_sa(input_DEM, clipmask_shp, output_DEM_ite)

        pixels = row.getValue("pixels")
        
        if Input_scaling_scheme == "Stone":
            # Process calculations on pixels
            ppg, pressure, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm \
                = processNumpy(output_DEM_ite)[:-1]  # No mask, we don't need the last value returned
            # Calculation of average production rates in the drainage basins:
            row.setValue("LAT_ave", SUM_LAT / ppg)
            row.setValue("LONG_ave", SUM_LONG / ppg)
            row.setValue("Alt_ave", SUM_DEM / ppg)
            row.setValue("Pneu_"+ input_correction_type, SUM_Pn / pixels)
            row.setValue("Psm_"+ input_correction_type, SUM_Psm / pixels)
            row.setValue("Pfm_"+ input_correction_type, SUM_Pfm / pixels)
            row.setValue("SFn_" + input_correction_type, SUM_SFn / pixels)
            row.setValue("SFf_" + input_correction_type, SUM_SFfm / pixels)
            row.setValue("SFs_" + input_correction_type, SUM_SFsm / pixels)
            row.setValue("pixfreeI",  ppg)
            row.setValue("Pres", pressure / ppg)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)

        elif Input_scaling_scheme == "LSD":
            # Process calculations on pixels
            if If_LSD_polynom_check:
                arcpy.AddMessage("LSD: Polynomial simplification")
                ppg, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm = processNumpy_LSD_poly(output_DEM_ite)[:-1]

            else:
                ppg, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn1, SUM_Psm1, SUM_Pfm1, SUM_SFn1, SUM_SFsm1, SUM_SFfm1 = processNumpy_LSD(output_DEM_ite)[:-1]
                SUM_Pn=SUM_Pn1[0]
                SUM_Psm=SUM_Psm1[0]
                SUM_Pfm=SUM_Pfm1[0]
                SUM_SFn=SUM_SFn1[0]
                SUM_SFsm=SUM_SFsm1[0]
                SUM_SFfm=SUM_SFfm1[0]
                
            row.setValue("LAT_ave", SUM_LAT / ppg)
            row.setValue("LONG_ave", SUM_LONG / ppg)
            row.setValue("Alt_ave", SUM_DEM / ppg)
            row.setValue("Pneu_"+ input_correction_type, SUM_Pn / pixels)
            row.setValue("Psm_"+ input_correction_type, SUM_Psm / pixels)
            row.setValue("Pfm_"+ input_correction_type, SUM_Pfm / pixels)
            row.setValue("SFn_" + input_correction_type, SUM_SFn / pixels)
            row.setValue("SFf_" + input_correction_type, SUM_SFfm / pixels)
            row.setValue("SFs_" + input_correction_type, SUM_SFsm / pixels)
            row.setValue("Nuclide", codenucl)
            row.setValue("Scaling", codescal)
            row.setValue("Options", optioncode)
            row.setValue("pixfreeI",  ppg)

        rows.updateRow(row)

    return output_shp

# Main

# =============================================================================================================
# Process: If the geological mask is checked, create new polygon shape file from the initial one using the mask
# in order to extract only quartz zone polygons of each drainage basins.
# =============================================================================================================

correction_type = ""
if If_geol_mask_checked:
    correction_type += "G"
    arcpy.AddMessage("Geological mask")

    # =================================================================================
    # Process : If the topographic mask is checked, in addition to the geological mask.
    # =================================================================================

    if If_topo_mask_checked:
        correction_type += "T"
        arcpy.AddMessage("Topographic mask")

        Output_geolquartz_polygons = geol_topo(Input_shape_file_of_the_drainage_basins, Input_DEM_raster, Input_geol_shape_file, For_quartz_zone_extraction, Input_topo_Mask_raster, Output_raster_DEM_of_the_drainage_basins,
                                               Output_geolquartz_polygons2,correction_type,Input_nuclide,Input_scaling_scheme)

        if If_VDM_checked:
            processVDM(Output_geolquartz_polygons, VDM, correction_type)

        # ===================================================================================================================================
        # Process: If, in addition to the geological and topographic mask, the glims mask is also checked: create new polygon shape file from
        # the geological one previously obtain using the glims mask
        # in order to extract only ice free quartz zone polygons of each drainage basins.
        # ===============================================================================

        if If_glims_mask_checked:
            correction_type += "I"
            arcpy.AddMessage("Ice mask")

            Output_geolquartz_polygons =pixel_count(Output_geolquartz_polygons, Input_DEM_raster,Output_raster_DEM_of_the_drainage_basins)
            Output_glims_polygons = glims_topo(Output_geolquartz_polygons, Input_DEM_raster, Input_topo_Mask_raster, Input_glims_shape_file, Output_raster_DEM_of_the_drainage_basins, Output_glims_polygons2, correction_type,Input_nuclide,Input_scaling_scheme)

            if If_VDM_checked:
                processVDM(Output_glims_polygons, VDM, correction_type)
                #removeFields(Output_glims_polygons, ["P_GT", "Topo_fc", "SFn_GT", "SFf_GT", "SFs_GT", "P_GTP", "SFn_GTP", "SFf_GTP", "SFs_GTP", "pixfreeI"])

    # ==========================================================================
    # Process : If the topographic mask is NOT checked (only the geological one)
    # ==========================================================================

    else:
        arcpy.AddMessage("No Topographic mask")

        Output_geolquartz_polygons = geol_notopo(Input_shape_file_of_the_drainage_basins, Input_DEM_raster, Input_geol_shape_file, For_quartz_zone_extraction, Output_raster_DEM_of_the_drainage_basins, Output_geolquartz_polygons2,correction_type,Input_nuclide,Input_scaling_scheme)

        if If_VDM_checked:
            processVDM(Output_geolquartz_polygons, VDM, correction_type)

        # ================================================================================
        # Process: If in addition to the geological mask, the glims mask is also checked :
        # create new polygon shape file from the geological one previously obtain using the glims mask in
        # order to extract only ice free quartz zone polygons of each drainage basins.
        # ============================================================================

        if If_glims_mask_checked:
            correction_type += "I"
            arcpy.AddMessage("Ice mask")

            Output_geolquartz_polygons =pixel_count(Output_geolquartz_polygons, Input_DEM_raster,Output_raster_DEM_of_the_drainage_basins)
            Output_glims_polygons = glims_notopo(Output_geolquartz_polygons, Input_DEM_raster, Input_glims_shape_file, Output_raster_DEM_of_the_drainage_basins, Output_glims_polygons2, correction_type,Input_nuclide,Input_scaling_scheme)

            if If_VDM_checked:
                processVDM(Output_glims_polygons, VDM, correction_type)
                #removeFields(Output_glims_polygons, ["P_G", "Topo_fc", "SFn_G", "SFf_G", "SFs_G", "P_GP", "SFn_GP", "SFf_GP", "SFs_GP", "pixfreeI"])

# ===============================================
# Process : If the geological mask is NOT checked
# ===============================================

else:
    arcpy.AddMessage("No Geological mask")

    # ========================================================
    # Process: If the topographic mask is the ONLY one checked
    # ========================================================

    if If_topo_mask_checked:
        correction_type += "T"
        arcpy.AddMessage("Topographic mask")
        # ========================================================================================================================
        # Process: Iterate the selection of all individual feature (i.e drainage basins) of the shape file of the drainage basins.
        # The goal is to extract as ascii file the DEM and topographic mask of each individual basin
        # in order to further calculate the production rates.
        # ===================================================

        Input_shape_file_of_the_drainage_basins = nogeol_topo(Input_shape_file_of_the_drainage_basins, Input_DEM_raster, Input_topo_Mask_raster, Output_raster_DEM_of_the_drainage_basins,correction_type,Input_nuclide,Input_scaling_scheme)

        if If_VDM_checked:
            processVDM(Input_shape_file_of_the_drainage_basins, VDM, correction_type)

        # =================================================================================
        # Process: If in addition to the topographic mask, the glims mask is also checked :
        # create new polygon shape file from the Input_shape_file_of_the_drainage_basins
        # using the glims mask in order to extract only ice free zone polygons of each drainage basins.
        # =============================================================================================

        if If_glims_mask_checked:
            correction_type += "I"
            arcpy.AddMessage("Ice mask")

            Input_shape_file_of_the_drainage_basins =pixel_count(Input_shape_file_of_the_drainage_basins, Input_DEM_raster,Output_raster_DEM_of_the_drainage_basins)
            Output_glims_polygons = glims_topo(Input_shape_file_of_the_drainage_basins, Input_DEM_raster, Input_topo_Mask_raster, Input_glims_shape_file, Output_raster_DEM_of_the_drainage_basins, Output_glims_polygons2, correction_type,Input_nuclide,Input_scaling_scheme)

            if If_VDM_checked:
                processVDM(Output_glims_polygons, VDM, correction_type)
                #removeFields(Output_glims_polygons, ["P_T", "Topo_fc", "SFn_T", "SFf_T", "SFs_T", "P_TP", "SFn_TP", "SFf_TP", "SFs_TP", "pixfreeI"])

    # ==========================
    # Process: If NO mask at all
    # ==========================

    else:
        arcpy.AddMessage("No Topographic mask")
        # ======================================================================================================
        # Process: Iterate the selection of all individual feature (i.e drainage basins) of the drainage basins.
        # The goal is to extract as ascii file the DEM of each individual basin in order to further calculate the production rates.
        # =========================================================================================================================

        Input_shape_file_of_the_drainage_basins = nogeol_notopo(Input_shape_file_of_the_drainage_basins, Input_DEM_raster, Output_raster_DEM_of_the_drainage_basins,Input_nuclide,Input_scaling_scheme)

        if If_VDM_checked:
            processVDM(Input_shape_file_of_the_drainage_basins, VDM, correction_type)

        # ==================================================================================================================================
        # Process: If the glims mask is the ONLY one checked: create new polygon shape file from the Input_shape_file_of_the_drainage_basins
        # using the glims mask in order to extract only ice free zone polygons of each drainage basins.
        # =============================================================================================

        if If_glims_mask_checked:
            correction_type += "I"
            arcpy.AddMessage("Ice mask")

            Input_shape_file_of_the_drainage_basins =pixel_count(Input_shape_file_of_the_drainage_basins, Input_DEM_raster,Output_raster_DEM_of_the_drainage_basins)
            Output_glims_polygons = glims_notopo(Input_shape_file_of_the_drainage_basins, Input_DEM_raster, Input_glims_shape_file, Output_raster_DEM_of_the_drainage_basins, Output_glims_polygons2, correction_type,Input_nuclide,Input_scaling_scheme)

            if If_VDM_checked:
                processVDM(Output_glims_polygons, VDM, correction_type)

# Deletes the temp_folder and all the associated contents
if os.path.exists(cur_dir + "/temp_folder"):
    shutil.rmtree(cur_dir + "/temp_folder")

# If the script is executed while the attribute table is open, fields are not automatically updated.
# There is no working arcpy function to reload the fields, this is a trick to make the refresh automatic.
# A new field is created and immediately removed to force this refresh.
arcpy.AddField_management(Input_shape_file_of_the_drainage_basins, "ReloadRows", "TEXT")
arcpy.DeleteField_management(Input_shape_file_of_the_drainage_basins, "ReloadRows")
