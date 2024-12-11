
# -*- coding: utf-8 -*-
#===============================================================================================================================================================================================
#BASINGA for QGIS: script which calculates the average cosmogenic production rate for individual drainage basin
#================================================================================================================================================================================================

#===========================================================================================================================================================================================
#PARMATERS AND PYTHON LIBRAIRIES
#===========================================================================================================================================================================================
#PYTHON LIBRARIES
#Open of the different python libraries needed later in the program:

# Import from plugin
from .parameters_definition import *
from .lsd import *

# Import from qgis core and existing functions
from qgis.core import *
from PyQt5.QtCore import *
import processing
import qgis.utils
from osgeo import gdal

# Tools for scientific calculation and files manipulation
import numpy as np
import math
from math import exp
from math import log
from scipy import interpolate
pi=math.pi

# Os tools
#from itertools import izip
import os
import os.path
from os import chdir, getcwd
import glob

VDM="MUSCHELER"

#===========================================================================================================================================================================================
#Auxiliary functions needed for the calculations
#===========================================================================================================================================================================================

def list_fields(type_):
    # Return a list of the names of the fields used in the process, requires a string specifying the type of corrections used : "GTI", "TI", "T"...
    # G stands for geological correction
    # T for Topographic shielding correction
    # I for Ice cover

    # Each field name is followed by "Double" which specify the kind of value that can be stored in the field.
    if type_ == 'GTI':
        list_new_fields = [["Pn_GTI","Double"],["SFn_GTI","Double"], ["Psm_GTI","Double"], ["SFsm_GTI","Double"], ["Pfm_GTI","Double"], ["SFfm_GTI","Double"],["topo_fc","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["LON_ave","Double"],["ALT_ave","Double"], ["pixels","Double"],["scaling","Double"],["nuclide","Double"],["options","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type_ == 'GI':
        list_new_fields = [["Pn_GI","Double"],["SFn_GI","Double"], ["Psm_GI","Double"], ["SFsm_GI","Double"], ["Pfm_GI","Double"], ["SFfm_GI","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"], ["LON_ave","Double"],["pixels","Double"],["scaling","Double"],["nuclide","Double"],["options","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type_ == 'TI':
        list_new_fields = [["Pn_TI","Double"],["SFn_TI","Double"], ["Psm_TI","Double"], ["SFsm_TI","Double"], ["Pfm_TI","Double"], ["SFfm_TI","Double"],["topo_fc","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"],["LON_ave","Double"], ["pixels","Double"],["scaling","Double"],["nuclide","Double"],["options","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type_ == 'I':
        list_new_fields = [["Pn_I","Double"],["SFn_I","Double"], ["Psm_I","Double"], ["SFsm_I","Double"], ["Pfm_I","Double"], ["SFfm_I","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"], ["LON_ave","Double"],["pixels","Double"],["scaling","Double"],["nuclide","Double"],["options","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type_ == 'GT':
        list_new_fields = [["Pn_GT","Double"],["SFn_GT","Double"], ["Psm_GT","Double"], ["SFsm_GT","Double"], ["Pfm_GT","Double"], ["SFfm_GT","Double"],["topo_fc","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"],["LON_ave","Double"], ["pixels","Double"],["scaling","Double"],["nuclide","Double"],["options","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type_ == 'G':
        list_new_fields = [["Pn_G","Double"],["SFn_G","Double"], ["Psm_G","Double"], ["SFsm_G","Double"], ["Pfm_G","Double"], ["SFfm_G","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"], ["LON_ave","Double"],["pixels","Double"],["scaling","Double"],["nuclide","Double"],["options","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type_ == 'T':
        list_new_fields = [["Pn_T","Double"],["SFn_T","Double"], ["Psm_T","Double"], ["SFsm_T","Double"], ["Pfm_T","Double"], ["SFfm_T","Double"],["topo_fc","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"],["LON_ave","Double"], ["pixels","Double"],["scaling","Double"],["nuclide","Double"],["options","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    else :
        list_new_fields = [["Pn_","Double"],["SFn_","Double"], ["Psm_","Double"], ["SFsm_","Double"], ["Pfm_","Double"], ["SFfm_","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"],["LON_ave","Double"], ["pixels","Double"],["scaling","Double"],["nuclide","Double"],["options","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]
    
    return list_new_fields

def code_type(type_):
    if type_ == 'GTI':
        code_type=1
    elif type_ == 'GI':
        code_type=2
    elif type_ == 'TI':
        code_type=3
    elif type_ == 'I':
        code_type=4
    elif type_ == 'GT':
        code_type=5
    elif type_ == 'G':
        code_type=6
    elif type_ == 'T':
        code_type=7
    else :
        code_type=8
    return code_type

def checkfield(inSHP,champs):
    #check wether or not the inSHP file already includes the fields needed to store the results of the requested calculations  as function of the selected options
    # list_field_u2 = [field.name() for field in inSHP.pendingFields()]
    list_field_u2 = [field.name() for field in inSHP.fields()] #MZ
    list_field02=[x.encode('UTF8') for x in list_field_u2]
    result=0
    if champs in list_field02:
        result=1      
    return result
    
def addFields(inSHP,listFields):
    # Add new field to a vector layer
    # inSHP must be an open layer, listFields a list of the new field names and their type : Int or Double
    # listFields = [[name,type],["newfield","Int"]]
    caps = inSHP.dataProvider().capabilities()
    if caps & QgsVectorDataProvider.AddFeatures:
        for field in listFields: # add each field specified, the type is respected
            check_field=checkfield(inSHP,field[0])
            if check_field ==0:
                if field[1] == "Int":
                    inSHP.dataProvider().addAttributes([QgsField(field[0], QVariant.Int)])
                elif field[1] == "Double":
                    inSHP.dataProvider().addAttributes([QgsField(field[0], QVariant.Double)])                                                                                                                                      

def assignSplitNumber (inSHP,field="split",type_ = "Int"):
    # assign a split number by creating a split flied and counting the polygons/features, inSHP must be an open layer, a field name and a value type can be specified

    addFields(inSHP,[[field, type_]])

    iter = inSHP.getFeatures() # get a list of the feature/polygon objects
    i = 1
    inSHP.startEditing() # open the editor mode
    # idx = inSHP.fieldNameIndex(field) # get the index of the field "field"
    idx = inSHP.fields().indexOf(field) #MZ
    for feature in iter:
        inSHP.changeAttributeValue(feature.id(), idx, i) # set i as value for the field idx of the polygon feature
        i += 1
    inSHP.commitChanges() # save the changes

def process_expression (expression):
    # extract variable (field), equ symbole (= , <, ...) and value from the expression and returns them
    express = expression.split('\"')
    var = express[1]
    express = express[2].split('\'')
    equ = express[0]
    val = express[1]
    # list_symbols = ['=', '!=', '>', '>=', '<', '<=', 'START_With', 'END_With']
    # Those are the symbols that can used later for the extraction
    if equ == '=':
        equ = 0
    elif equ == '!=':
        equ = 1
    elif equ == '>':
        equ = 2
    elif equ == '>=':
        equ = 3
    elif equ == '<':
        equ = 4
    elif equ == '<=':
        equ = 5
    elif equ == 'START_With':
        equ = 6
    elif equ == 'END_With':
        equ = 7

    return var, equ, val

def openLayer(type_,path,legendName,pilot='ogr'):
    #open vector or raster layer, requires the type ("raster" or "layer"), the path of the wanted file, the name that will be used
    # in the qgis legend and the pilot for vector file (defautl set as "ogr", no change needed)
    # Return the opened layer as an qgis layer object
    if type_ == 'raster':
        objectName = QgsRasterLayer(path, legendName)
    elif type_ == 'vector':
        objectName = QgsVectorLayer(path, legendName, pilot)

    QgsMapLayerRegistry.instance().addMapLayer(objectName) #Add the loaded layer in qgis interface

    return objectName

def getExtent(inSHP):
    # return extent of an open layer as a string "x,x,x,x", requires an open layer
    e = inSHP.extent().toString() # get the extent and transform it to a string value
    e = e.replace(' : ', ',')
    e = e.encode('ascii', 'ignore')
    e = e.split(',')
    ext = e[0] + ',' + e[2] + ',' + e[1] + ',' + e[3]

    return ext

def createPolyNames(shp_file):
    # from a file path/name return a list of path/name which correspond to the same name but with different extension

    names = []

    names.append(shp_file[:-4] + '.tif')
    names.append(shp_file[:-4] + 'mask' + '.tif')
    names.append(shp_file[:-4] + '.asc')
    names.append(shp_file[:-4] + '_ts.asc')

    return names

def selectSaveByAttribute(inSHP,expression,outPath) :
    # select by attribute and save the selection, requires an open vector layer , the
    # expression to make the selection and a path for the outcome
    var, equ, val = process_expression(expression)

    processing.runalg('qgis:selectbyattribute', inSHP, var, equ, val) #qgis command to selectbyAttribute

    processing.runalg('qgis:saveselectedfeatures', inSHP, outPath) #qgis command to save the selection

def splitSHP(inSHP,fieldName,dir):
    # split an open layer in dir according to a field name and return the list of the paths of the created shp files
    # processing.runalg('qgis:splitvectorlayer', inSHP, fieldName, dir) # qgis command to split a vector layer to x files containing each one the polygons sharing a common field value
    #MZ
    processing.run("native:splitvectorlayer", {
        'INPUT': inSHP,
        'FIELD': fieldName,
        'OUTPUT': dir
    })
    # Make a list of every shp file in the temp directory that was created above (it avoids any other shapefile file in the folder)
    temp_dir_access = dir + '/' + '*.shp'
    list_shp_global = glob.glob(temp_dir_access) #
    list_shp = []
    for shpfile in list_shp_global:
        if shpfile.find(fieldName) != -1:
            list_shp.append(shpfile)
    return list_shp

def clipRaster_byLayer(inDEM,layer,outDEM):
    # clip a raster layer (opened) by a vector layer and save it in the file outDEM
    # It uses a GDAL script which is different in qgis 2.14 and 2.8

    # We get the qgis version to pick the right command line
    # v = qgis.utils.QGis.QGIS_VERSION
    v = Qgis.QGIS_VERSION ##MZ
    v = v.split(".")

    # Then we perform the gdal script
    if v[0] == "2" and v[1] == "8":
        processing.runalg('gdalogr:cliprasterbymasklayer', inDEM, layer, None, False, False, None, outDEM)
    elif v[0] =="2" and v[1] == "14":
        processing.runalg('gdalogr:cliprasterbymasklayer', inDEM, layer, None, False, True, True, 5, 4, 75, 6, 1, False, 0, False, None, outDEM)
    elif v[0] == "2" and v[1] == "16":
        processing.runalg('gdalogr:cliprasterbymasklayer', inDEM, layer, None, False, True, True, 5, 4, 75, 6, 1, False, 0, False, None, outDEM)
    else:
        print("the plugin is may no be able to deal with your version of qgis : the gdal tool cliprasterbymasklayer is different")
        
    return outDEM

def readHeaderfile(file):
    # read the header of one opened file and return the selected information
    
    # How do read the ascii files?#
    # ncols         3401
    # nrows         3264
    # xllcorner     82.828911273175
    # yllcorner     41.747419874025
    # dx            0.000277777778 or cellsize = ....
    # dy            0.000277690002
    # NODATA_value  -9999
    # -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999
    # -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999
    # -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999 -9999

    # Read the header and its informations:
    # read the number of columns in line 1
    line = file.readline()
    title, ncols = line.split()
    ncols = int(ncols)
    # read the number of lines in line 2
    line = file.readline()
    title, nrows = line.split()
    nrows = int(nrows)
    # read the longitude coordinates of the left lower corner in line 3
    line = file.readline()
    title, xllcorner = line.split()
    xllcorner = float(xllcorner)
    # read the latitude coordinates of the left lower corner in line 4
    line = file.readline()
    title, yllcorner = line.split()
    yllcorner = float(yllcorner)

    # read the cellsize / dx of the DEM
    # make sure that either cellsize or dy is saved
    line = file.readline()
    title, cellsize = line.split()
    cellsize = float(cellsize)

    if title == 'cellsize':
        dy = cellsize

    else:
        # read the dy cellsize of the DEM
        line = file.readline()
        title, dy = line.split()
        dy = float(dy)
    # read the non data values
    line = file.readline()
    title, NODATA_value = line.split()
    
    return ncols, nrows, xllcorner, yllcorner, dy, NODATA_value

def getFieldIndex(inSHP,listFieldName):
    # return a list of field index from a opened shp file and a list of field names
    list_fields_index = []
    for field in listFieldName:
        # idx = inSHP.fieldNameIndex(field)
        idx = inSHP.fields().indexOf(field) #MZ
        list_fields_index.append(idx)

    return list_fields_index


def updateResults(inSHP, listFieldName, results, split='split'):
    """
    MZ
    Updates the attribute table of a layer based on the results matrix.

    Args:
        inSHP (QgsVectorLayer): The input shapefile layer.
        listFieldName (list): List of field names to update.
        results (list): A matrix of results to copy to the layer.
        split (str): The field name used to match features (default: 'split').

    Returns:
        None
    """
    list_fields_index = getFieldIndex(inSHP, listFieldName)

    if inSHP.fields().indexOf(split) == -1:
        print(f"Error: Field '{split}' not found in the layer.")
        return

    iter_features = inSHP.getFeatures()
    inSHP.startEditing()

    for feature in iter_features:
        num = feature[split]
        t = -1
        it = 0

        while t == -1:
            if it >= len(results):
                print(f"Warning: No match found for feature with '{split}' = {num}. Skipping.")
                break

            if len(results[it]) > 0:  # Ensure results[it] exists and is not empty
                num2 = results[it][0]
                if num == num2:
                    t = it
            it += 1

        if t == -1:
            continue

        for j in range(len(list_fields_index)):
            if t < len(results) and j + 1 < len(results[t]):
                inSHP.changeAttributeValue(feature.id(), list_fields_index[j], float(results[t][j + 1]))
            else:
                print(f"Warning: Missing data for feature '{split}' = {num}'. Skipping field update.")

    inSHP.commitChanges()
    print(f"Attribute table updated successfully.")

# def updateResults(inSHP,listFieldName,results,split='split'):
#     # requires an opened layer, the list of the field names used, the matrix of the results, and the name of the field used to count the polygons
#     list_fields_index = getFieldIndex(inSHP, listFieldName) # get the list the list of the field indexes to work with it
#
#     iter = inSHP.getFeatures()
#     inSHP.startEditing()
#
#     # copy the values in results in the matching place in the attribute table of inSHP
#     i = 0
#     for features in iter:
#         num = features[split]
#         t = -1
#         it = 0
#         while t == -1:
#             num2 = results[it][0]
#             if num == num2: # Make sure that the split number match to that the two table present the same line when copying
#                 t = it
#             it += 1
#         for j in range(len(list_fields_index)):
#             inSHP.changeAttributeValue(features.id(), list_fields_index[j], float(results[t][j + 1])) # we save the results in the wanted fields by specifying the feature, the index of the field and the value
#         i += 1
#
#     inSHP.commitChanges()

def deleteFields(inSHP,listFields):
    # requires an opened layer and a list of names of the fields to delete
    list_fields_index_todelete = getFieldIndex(inSHP,listFields) # get their indexes
    caps = inSHP.dataProvider().capabilities()
    if caps & QgsVectorDataProvider.DeleteAttributes:
        inSHP.dataProvider().deleteAttributes(list_fields_index_todelete) # delete accordingly to their indexes

    inSHP.updateFields()

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
    
    non_zero = latitude > 0
    lat_min=np.min(latitude[non_zero])
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

def pixel_count(inSHP,inDEM):

    polygonNames = createPolyNames(inSHP)
    clipRaster_byLayer(inDEM, inSHP, polygonNames[0])

    DEM_loaded = openLayer('raster',polygonNames[0], inSHP[:-4] + 'dtm')
    extent = getExtent(DEM_loaded)
    QgsMapLayerRegistry.instance().removeMapLayer(DEM_loaded.id())

    processing.runalg("gdalogr:translate", polygonNames[0], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[2])

    f_DEM = open(polygonNames[2], 'r') # open the file in reading only mode
    DEM_array=[]    
    ncols, nrows, left_long, bot_lat, cellsize, NODATA_value = readHeaderfile(f_DEM) # read and extract metadata from the header...to skip the first lines
    # Reading data and transform to a numpy array
    for line_DEM in f_DEM:
        # splitting the line according to the tabs
        list_line_DEM = line_DEM.split()
        # transforming each element of the list into numbers
        list_DEM = [(float(x)) for x in list_line_DEM]
        DEM_array_line=np.asarray(list_DEM)
        DEM_array.append(DEM_array_line)
    DEM_array=np.asarray(DEM_array)
    pixel_with_data = DEM_array >= 0
    pixels= np.sum(DEM_array > 0)

    f_DEM.close()

    return pixels

#===========================================================================================================================================================================================
#Function to calculate the scaling factors and production rates accordingly either to the Lal/St or to the LSD model
#===========================================================================================================================================================================================
#LAL/STONE:
def processNumpy(DEM,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,MASK=None):
    """
    Computes, according to the Lal/Stone model the sum of the scaling factors, production rates and other parameters
    Inputs:
        DEM            : raster file
        MASK(optional) : topographic mask
        SLHL for the different produciton pathways
    Outputs:
        SUM_Ben  : neutrons production rates - g/at/yr
        SUM_SFn  : neutrons Stone factors
        SUM_Besm : slow muons production rates - g/at/yr
        SUM_SFsm : slow muons Stone factors
        SUM_Befm : fast muons production rates - g/at/yr
        SUM_SFfm : fast muons Stone factors
        SUM_P    : pressures - hPa
        SUM_LAT  : latitudes - degrees
        SUM_DEM  : altitudes of the DEM - m
        pixels   : number of pixels with data
        SUM_MASK : topographic factors (if selected)
    Note:
        - All the outputs are sums over all pixels.
        - If a pixel from the DEM is negative, it will be seen as a NoData value and will not be counted. This may lead to unexpected errors.
        - The data pixels from the mask has to fit exactly with the data pixels from the DEM.
    """
    f_DEM = open(DEM, 'r') # open the file in reading only mode
    if MASK is not None:
        f_MASK = open(MASK, 'r')
        ncols2, nrows2, left_long2, bot_lat2, cellsize2, NODATA_value2 = readHeaderfile(f_MASK) # read and extract metadata from the header of the MASK (not needed but just to pass the header lines)
    
    ncols, nrows, left_long, bot_lat, cellsize, NODATA_value = readHeaderfile(f_DEM) # read and extract metadata from the header
    
    top_lat = bot_lat + (nrows * cellsize)
    right_long=left_long+(ncols*cellsize)
    top_lat = float(top_lat)

    DEM_array=[]
    MASK_array=[]

    # Reading data and transform to a numpy array
    for line_DEM in f_DEM:

        # splitting the line according to the tabs
        list_line_DEM = line_DEM.split()

        # transforming each element of the list into numbers
        list_DEM = [(float(x)) for x in list_line_DEM]
        DEM_array_line=np.asarray(list_DEM)
        DEM_array.append(DEM_array_line)

    DEM_array=np.asarray(DEM_array)

    pixel_with_data = DEM_array >= 0  # Boolean array that stores where there is Data pixel
    
    taille0=np.shape(DEM_array)

    if MASK is not None:  # Test if there is a mask
        for line_MASK in f_MASK:

        # splitting the line according to the tabs
            list_line_MASK = line_MASK.split()

        # transforming each element of the list into numbers
            list_MASK = [(float(x)) for x in list_line_MASK]
            MASK_array_line=np.asarray(list_MASK)
            MASK_array.append(MASK_array_line)
        MASK_array=np.asarray(MASK_array)
        
    else:  # There is no MASK, MASK_array is the same shape than the DEM_array, filled with ones
        MASK_array = np.ones(DEM_array.shape)

    if taille0[0] != nrows:#correction made for bug in the function RasterToNumpy which sometime may generate additional rows when compared to the true original extracted DEM
        nrows=taille0[0]

    # Sea level pressure, sea level temperature and adiabatic lapse rate at the center of the ASCII file from the ERA40 atmosphere model
    center_LAT = bot_lat + 0.5 * nrows * cellsize
    center_LONG = left_long + 0.5 * ncols * cellsize
    sea_level_P, sea_level_T, dtdz = ERA40_Data(center_LAT, center_LONG)

    a, b, c, d, e = DEMtoStoneCoefficients(DEM_array, bot_lat, cellsize)

    LAT = np.linspace(bot_lat, top_lat, nrows).reshape(nrows, 1) * pixel_with_data  # Latitudes array with zeros where DEM is negative
    
    rate=correction_latitude(LAT,cellsize)#weighting for latitudinal variation of pixel surface

    # PRESSURE TABLE as a function of elevation (Stone, 2000):
    # Nodata_values are excluded
    P = sea_level_P * np.exp(-0.03417 / dtdz * (np.log(sea_level_T) - np.log(sea_level_T - dtdz * DEM_array[pixel_with_data])))

    SF_Be_neutron = (a + b * np.exp(-P / 150) + c * P + d * P ** 2 + e * P ** 3) * MASK_array[pixel_with_data]
    ProdDEM_Be_neutron = SLHLP_Be_neutron * SF_Be_neutron
    SF_Be_neutron=SF_Be_neutron*rate
    ProdDEM_Be_neutron=ProdDEM_Be_neutron*rate

    SF_Be_fastmuon = np.exp((1013.25 - P) / at_attenuation_fastmuon) * MASK_array[pixel_with_data]
    ProdDEM_Be_fastmuon = SLHLP_Be_fastmuon * SF_Be_fastmuon
    SF_Be_fastmuon = SF_Be_fastmuon *rate
    ProdDEM_Be_fastmuon=ProdDEM_Be_fastmuon*rate

    SF_Be_slowmuon = np.exp((1013.25 - P) / at_attenuation_slowmuon) * MASK_array[pixel_with_data]
    ProdDEM_Be_slowmuon = SLHLP_Be_slowmuon * SF_Be_slowmuon
    SF_Be_slowmuon = SF_Be_slowmuon *rate
    ProdDEM_Be_slowmuon=ProdDEM_Be_slowmuon*rate

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
        orint("<!> No pixel with data found in the basin <!> There must a polygon with zero pixel size. Please check your basins and the number of polygons")

    f_DEM.close()
    if MASK is not None:
        f_MASK.close() # closing the files

    return pixels, SUM_Pre, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK

#LSD from polynomial interpolation of Lal/Stone factors
def processNumpy_LSD_poly(DEM,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide,MASK=None):
    """
    Computes, using a polynomial interpolation form Lal/Stone factors, the sum of the scaling factors, production rates according to the LSD model
    Inputs:
        DEM            : raster file
        MASK(optional) : topographic mask
        SLHL for the different produciton pathways
    Outputs:
        SUM_Ben  : neutrons production rates - g/at/yr
        SUM_SFn  : neutrons Stone factors
        SUM_Besm : slow muons production rates - g/at/yr
        SUM_SFsm : slow muons Stone factors
        SUM_Befm : fast muons production rates - g/at/yr
        SUM_SFfm : fast muons Stone factors
        SUM_P    : pressures - hPa
        SUM_LAT  : latitudes - degrees
        SUM_DEM  : altitudes of the DEM - m
        pixels   : number of pixels with data
        SUM_MASK : topographic factors (if selected)
    Note:
        - All the outputs are sums over all pixels.
        - If a pixel from the DEM is negative, it will be seen as a NoData value and will not be counted. This may lead to unexpected errors.
        - The data pixels from the mask has to fit exactly with the data pixels from the DEM.
    """
    import random
    f_DEM = open(DEM, 'r') # open the file in reading only mode
    if MASK is not None:
        f_MASK = open(MASK, 'r')
        ncols2, nrows2, left_long2, bot_lat2, cellsize2, NODATA_value2 = readHeaderfile(f_MASK) # read and extract metadata from the header of the MASK (not needed but just to pass the header lines)
    
    ncols, nrows, left_long, bot_lat, cellsize, NODATA_value = readHeaderfile(f_DEM) # read and extract metadata from the header
    
    top_lat = bot_lat + (nrows * cellsize)
    right_long=left_long+(ncols*cellsize)
    top_lat = float(top_lat)

    DEM_array=[]
    MASK_array=[]

    # Reading data and transform to a numpy array
    for line_DEM in f_DEM:

        # splitting the line according to the tabs
        list_line_DEM = line_DEM.split()

        # transforming each element of the list into numbers
        list_DEM = [(float(x)) for x in list_line_DEM]
        DEM_array_line=np.asarray(list_DEM)
        DEM_array.append(DEM_array_line)

    DEM_array=np.asarray(DEM_array)
  
    pixel_with_data = DEM_array >= 0  # Boolean array that stores where there is Data pixel
    
    DEM_bon=DEM_array[pixel_with_data]

    taille0=np.shape(DEM_array)

    if MASK is not None:  # Test if there is a mask
        for line_MASK in f_MASK:

        # splitting the line according to the tabs
            list_line_MASK = line_MASK.split()

        # transforming each element of the list into numbers
            list_MASK = [(float(x)) for x in list_line_MASK]
            MASK_array_line=np.asarray(list_MASK)
            MASK_array.append(MASK_array_line)
        MASK_array=np.asarray(MASK_array)
    else:  # There is no MASK, MASK_array is the same shape than the DEM_array, filled with ones
        MASK_array = np.ones(DEM_array.shape)

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
    
    SF_Be_neutron = (a + b * np.exp(-P / 150) + c * P + d * P ** 2 + e * P ** 3) 

    SF_Be_fastmuon = np.exp((1013.25 - P) / at_attenuation_fastmuon) 

    SF_Be_slowmuon = np.exp((1013.25 - P) / at_attenuation_slowmuon) 

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
    
    for o in range(1,10):
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

    for i, dem in enumerate(DEM_extract_end):  # Loop over the DEM array, because LSD code is not fully vectorized and needs single values as inputs
        if dem <= 0: continue # Skip DEM negative values
            
        tmp_time_vec, tmp_SF_vecSP, tmp_SF_vecMU= LSDv9(Lat_extract_end[i], Lon_extract_end[i], dem, 0, 0, -1, Input_nuclide, MUSCHELERm)
     
        lsd_sp_extract_end.append(tmp_SF_vecSP[0])
        lsd_mu_extract_end.append(tmp_SF_vecMU[0])

    ##Step 4: polynomial fit between Stone and LSD factors on the sampled pixels:

    Preg = np.polyfit(St_sp_extract_end, lsd_sp_extract_end, 4)# determine the polynom to calculte SP LSD production from stone production
    Preg2 = np.polyfit(St_sm_extract_end, lsd_mu_extract_end, 4)# determine the polynom to calculte slow muons LSD production from Braucher production
    Preg3 = np.polyfit(St_fm_extract_end, lsd_mu_extract_end, 4)# determine the polynom to calculte fast muons LSD production from Braucher production

    ##Step 5:extrapollation of  LSD factors to all cells based on the polifit relation:
    #Total_Prod_LSD = []
    #for i in range(len(Total_Prod_st)):
    #    Total_Prod_LSD.append(Preg[0] * Total_Prod_st[i] ** 4 + Preg[1] * Total_Prod_st[i] ** 3 + Preg[2] * Total_Prod_st[i] ** 2 + Preg[3] * Total_Prod_st[i] + Preg[4])

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

    f_DEM.close()
    if MASK is not None:
        f_MASK.close() # closing the files

    return pixels,SUM_LAT, SUM_LONG, SUM_DEM,SUM_Ben, SUM_Besm, SUM_Befm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK

#LSD
def processNumpy_LSD(DEM,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide,MASK=None):
    """
    Computes, according to the LSD model the sum of the scaling factors, production rates and other parameters
    Inputs:
        DEM            : raster file
        MASK(optional) : topographic mask
        SLHL for the different produciton pathways
    Outputs:
        SUM_Ben  : neutrons production rates - g/at/yr
        SUM_SFn  : neutrons Stone factors
        SUM_Besm : slow muons production rates - g/at/yr
        SUM_SFsm : slow muons Stone factors
        SUM_Befm : fast muons production rates - g/at/yr
        SUM_SFfm : fast muons Stone factors
        SUM_P    : pressures - hPa
        SUM_LAT  : latitudes - degrees
        SUM_DEM  : altitudes of the DEM - m
        pixels   : number of pixels with data
        SUM_MASK : topographic factors (if selected)
    Note:
        - All the outputs are sums over all pixels.
        - If a pixel from the DEM is negative, it will be seen as a NoData value and will not be counted. This may lead to unexpected errors.
        - The data pixels from the mask has to fit exactly with the data pixels from the DEM.
    """
    f_DEM = open(DEM, 'r') # open the file in reading only mode
    if MASK is not None:
        f_MASK = open(MASK, 'r')
        ncols2, nrows2, left_long2, bot_lat2, cellsize2, NODATA_value2 = readHeaderfile(f_MASK) # read and extract metadata from the header of the MASK (not needed but just to pass the header lines)
    
    ncols, nrows, left_long, bot_lat, cellsize, NODATA_value = readHeaderfile(f_DEM) # read and extract metadata from the header
    
    top_lat = bot_lat + (nrows * cellsize)
    right_long=left_long+(ncols*cellsize)
    top_lat = float(top_lat)

    DEM_array=[]
    MASK_array=[]

    # Reading data and transform to a numpy array
    for line_DEM in f_DEM:

        # splitting the line according to the tabs
        list_line_DEM = line_DEM.split()

        # transforming each element of the list into numbers
        list_DEM = [(float(x)) for x in list_line_DEM]
        DEM_array_line=np.asarray(list_DEM)
        DEM_array.append(DEM_array_line)

    DEM_array=np.asarray(DEM_array)

    pixel_with_data = DEM_array >= 0  # Boolean array that stores where there is Data pixel

    taille0=np.shape(DEM_array)

    if MASK is not None:  # Test if there is a mask
        
        for line_MASK in f_MASK:
        # splitting the line according to the tabs
            list_line_MASK = line_MASK.split()
        # transforming each element of the list into numbers
            list_MASK = [(float(x)) for x in list_line_MASK]
            MASK_array_line=np.asarray(list_MASK)
            MASK_array.append(MASK_array_line)
        MASK_array=np.asarray(MASK_array)
        
    else:  # There is no MASK, MASK_array is the same shape than the DEM_array, filled with ones
        MASK_array = np.ones(DEM_array.shape)

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
  
        tmp_time_vec, tmp_SF_vecSP, tmp_SF_vecMU = LSDv9(LAT.flatten()[i], LONG.flatten()[i], dem, 0, 0, -1, Input_nuclide, MUSCHELERm)
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

    f_DEM.close()
    if MASK is not None:
        f_MASK.close() # closing the files

    return pixels, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK

#=======================================================================================================================================================================================
#Functions for the paleomagnetic corrections
#=======================================================================================================================================================================================
def VDM_process_lsd(inSHP, ConcField,SLHLP_Be_neutron, SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide,type_):

    features = inSHP.getFeatures()

    # Create a list of the field used for calculation and saves
    if type_ == 'GTI':
        list_fields_in = ["Pn_GTI","SFn_GTI","LAT_ave", "LON_ave","ALT_ave", ConcField,"Time"]
    elif type_ == 'GI':
        list_fields_in = ["Pn_GI","SFn_GI","LAT_ave", "LON_ave","ALT_ave", ConcField,"Time"]
    elif type_ == 'TI':
        list_fields_in = ["Pn_TI","SFn_TI","LAT_ave", "LON_ave","ALT_ave", ConcField,"Time"]
    elif type_ == 'I':
        list_fields_in = ["Pn_I","SFn_I","LAT_ave", "LON_ave","ALT_ave", ConcField,"Time"]
    elif type_ == 'T':
        list_fields_in = ["Pn_T","SFn_T","LAT_ave","LON_ave", "ALT_ave", ConcField,"Time"]
    elif type_ == 'G':
        list_fields_in = ["Pn_G","SFn_G","LAT_ave","LON_ave", "ALT_ave", ConcField,"Time"]
    elif type_ == 'GT':
        list_fields_in = ["Pn_GT","SFn_GT","LAT_ave","LON_ave", "ALT_ave", ConcField,"Time"]
    else:
        list_fields_in = ["Pn_","SFn_","LAT_ave", "LON_ave","ALT_ave", ConcField,"Time"]
        # list index in [0    ,  1    ,  2     ,   3   ,      4      , 5           6]

    # get the list of field indexes
    list_index_in = []
    for field in list_fields_in:
        list_index_in.append(inSHP.fieldNameIndex(field))

    list_fields_out=["Pn_P","SFn_P","Psm_P","SFsm_P","Pfm_P","SFfm_P","Time"]
    #                  0       1       2       3        4       5       6
    list_index_out = []
    for field in list_fields_out:
        list_index_out.append(inSHP.fieldNameIndex(field))
    
    inSHP.startEditing()

    for feature in features:
        # we get the values for calculation, trying to pick the right ones depending on which correction is operated
        Pn = feature.attributes()[list_index_in[0]]
        Conc = feature.attributes()[list_index_in[5]]
        LAT = feature.attributes()[list_index_in[2]]
        LON = feature.attributes()[list_index_in[3]]
        Z = feature.attributes()[list_index_in[4]]

        Time = Conc / Pn

        # calculation
        Sneu_corr,Smu_corr=timecorrect_LSD(LAT,LON,Z,Time,nuclide)
        Pneu_corr=Sneu_corr*SLHLP_Be_neutron
        Psm_corr=Smu_corr*SLHLP_Be_slowmuon
        Pfm_corr=Smu_corr*SLHLP_Be_fastmuon

        # we create a list of field where we want to save the results
        list_values = [Pneu_corr,Sneu_corr,Psm_corr,Smu_corr,Pfm_corr,Smu_corr,Time]
        #                0           1         2       3           4       5    6

        # we save the results            
        inSHP.changeAttributeValue(feature.id(), list_index_out[0], float(list_values[0]))#Pn_P
        inSHP.changeAttributeValue(feature.id(), list_index_out[1], float(list_values[1]))#SFn_P
        inSHP.changeAttributeValue(feature.id(), list_index_out[2], float(list_values[2]))#Psm_P
        inSHP.changeAttributeValue(feature.id(), list_index_out[3], float(list_values[3]))#SFsm_P
        inSHP.changeAttributeValue(feature.id(), list_index_out[4], float(list_values[4]))#Pfm_P
        inSHP.changeAttributeValue(feature.id(), list_index_out[5], float(list_values[5]))#SFfm_P
        inSHP.changeAttributeValue(feature.id(), list_index_out[6], list_values[6])#Time

    inSHP.commitChanges()

def timecorrect_LSD(lat_ave,long_ave,Z_ave,AGE,nuclide):
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
    time_vec, scalingFactors_vecSP,scalingFactors_vecMU = LSDv9(lat_ave, long_ave, Z_ave, 0, AGE, -1, nuclide, MUSCHELERm)
    SFn_integre = LSDtimeInteg(10, time_vec, scalingFactors_vecSP)
    SFmu_integre = LSDtimeInteg(10, time_vec, scalingFactors_vecMU)
    
    return SFn_integre,SFmu_integre

def VDM_process(inSHP,VDM,ConcField,SLHLP_Be_neutron,type_):
    # from an opened layer, a specified VDM ("GLOPIS", "MUSCHELER" and "VDM2000"), the name of the field containing the concentration and the type of corretion previously applied)

    # import the correct VDM depending on the one choosen by the user
    if (VDM == "GLOPIS"):
        VDM_ = GLOPIS
    elif (VDM == "MUSCHELER"):
        VDM_ = MUSCHELER
    elif (VDM == "VDM2000"):
        VDM_ = VDM2000

    features = inSHP.getFeatures()

    # Create a list of the field used for calculation and record
    if type_ == 'GTI':
        list_fields_in = ["Pn_GTI","SFn_GTI","LAT_ave", "H_ave", ConcField,"Time"]
    elif type_ == 'GI':
        list_fields_in = ["Pn_GI","SFn_GI","LAT_ave", "H_ave", ConcField,"Time"]
    elif type_ == 'TI':
        list_fields_in = ["Pn_TI","SFn_TI","LAT_ave", "H_ave", ConcField,"Time"]
    elif type_ == 'I':
        list_fields_in = ["Pn_I","SFn_I","LAT_ave", "H_ave", ConcField,"Time"]
    elif type_ == 'T':
        list_fields_in = ["Pn_T","SFn_T","LAT_ave", "H_ave", ConcField,"Time"]
    elif type_ == 'G':
        list_fields_in = ["Pn_G","SFn_G","LAT_ave", "H_ave", ConcField,"Time"]
    elif type_ == 'GT':
        list_fields_in = ["Pn_GT","SFn_GT","LAT_ave", "H_ave", ConcField,"Time"]
    else:
        list_fields_in = ["Pn_","SFn_","LAT_ave", "H_ave", ConcField,"Time"]
        # list index in [0    ,  1    ,  2     ,   3   ,      4      , 5]

    # get the list of field indexes
    list_index_in = []
    for field in list_fields_in:
        list_index_in.append(inSHP.fieldNameIndex(field))

    list_fields_out=["Pn_P","SFn_P","Time"]
    list_index_out = []
    for field in list_fields_out:
        list_index_out.append(inSHP.fieldNameIndex(field))
    
    inSHP.startEditing()

    for feature in features:
        # we get the values for calculation, trying to pick the right ones depending on which correction is operated
        Pn = feature.attributes()[list_index_in[0]]
        Sn_non_corrected = feature.attributes()[list_index_in[1]]
        Conc = feature.attributes()[list_index_in[4]]
        LAT = feature.attributes()[list_index_in[2]]
        P = feature.attributes()[list_index_in[3]]

        Time = Conc / Pn

        # calculation
        Corr_fact = Rc_correctedStoneFactors(VDM_, LAT, P, Time)
        Sneu_corr=Sn_non_corrected*Corr_fact
        Pneu_corr=Sneu_corr*SLHLP_Be_neutron

        # we create a list of field where we want to save the results
        list_values = [Pneu_corr,Sneu_corr,Time]

        # we save the results            
        inSHP.changeAttributeValue(feature.id(), list_index_out[0], float(list_values[0]))
        inSHP.changeAttributeValue(feature.id(), list_index_out[1], float(list_values[1]))
        inSHP.changeAttributeValue(feature.id(), list_index_out[2], list_values[2])

    inSHP.commitChanges()

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

#============================================================================================================================================================================================================
#Main functions used to extract the D.E.M, the elevations, latitudes and calculate the production rates and scaling factors accordingly to the optionS selected (TS maks, geological mask, glacial coer, VDM)
#============================================================================================================================================================================================================
def glims_TS(shp, inDEM, temp_dir, in_ts,use_geol_mask, inGeolShp,qz_express,outQZshp, use_VDM, VDM, ConcField, inGLIMS, outGLIMS, types,lsORlsd,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide):
 
    if use_geol_mask == True:    
        temp = temp_dir + '/' + 'geolSHP_temp.shp'
        selectSaveByAttribute(inGeolShp, qz_express, temp)
        temp2 = temp_dir + '/' + 'geolSHP_temp2.shp'
        processing.runalg('qgis:clip', shp, temp, temp2)
        outSHP = openLayer('vector',temp2, temp2[:-4], "ogr")
        selectSaveByAttribute(outSHP,'\"pixels\">=\'0\'',outQZshp)
        inSHP=outQZshp
        inSHP = openLayer('vector', temp2, temp2[:-4], "ogr")
        assignSplitNumber(inSHP, "splitQz")
        list_shp_Qz=splitSHP(inSHP, "splitQz", temp_dir)
    else:
        inSHP=shp
        assignSplitNumber(inSHP, "split")
        list_shp_0=splitSHP(inSHP, "split", temp_dir)
   
    temp3 = temp_dir + '/' + 'out_ice_temp.shp'
    processing.runalg('qgis:difference',inSHP,inGLIMS, True,temp3)
    
    outSHP = openLayer('vector', temp3, temp3[:-4], "ogr")

    addFields(outSHP, list_fields(types))

    assignSplitNumber(outSHP, "splitIce")

    results = []

    # Create one file for each polygon of the shapefile
    list_shp = splitSHP(outSHP, "splitIce", temp_dir)
    length=len(list_shp)
    
    for i in range (0,length):
        shp_polygon=list_shp[i]
        polygonNames = createPolyNames(shp_polygon)

        clipRaster_byLayer(inDEM, shp_polygon, polygonNames[0])
        clipRaster_byLayer(in_ts, shp_polygon, polygonNames[1])


        DEM_loaded = openLayer('raster', polygonNames[0], shp_polygon[:-4] + 'dtm')
        extent = getExtent(DEM_loaded)
        QgsMapLayerRegistry.instance().removeMapLayer(DEM_loaded.id())

        processing.runalg("gdalogr:translate", polygonNames[0], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[2])
        processing.runalg("gdalogr:translate", polygonNames[1], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[3])
         
        if lsORlsd == "ls": 
            ppg, SUM_Pre, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK = processNumpy(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,polygonNames[3])        
        elif lsORlsd == "lsd":
            ppg, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide,polygonNames[3])
        else:
            ppg, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD_poly(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide,polygonNames[3])

        if use_geol_mask == True:
            shp_polygon_Qz=list_shp_Qz[i]
            pixels=pixel_count(shp_polygon_Qz,inDEM)
        else:
            shp_polygon_0=list_shp_0[i]
            pixels=pixel_count(shp_polygon_0,inDEM)
                    
        ProdDEM_Ben = SUM_Pn / pixels
        SF_Ben = SUM_SFn / pixels
        ProdDEM_Besm = SUM_Psm / pixels
        SF_Besm = SUM_SFsm / pixels
        ProdDEM_Befm = SUM_Pfm / pixels
        SF_Befm = SUM_SFfm / pixels
        Alt_average = SUM_DEM / pixels
        if lsORlsd == "ls":
            P_average=SUM_Pre/pixels
            scal=1
        else:
            P_average=0
            scal=0
        LAT_ave=SUM_LAT/pixels
        LON_ave=SUM_LONG/pixels
        Topo_fc=SUM_MASK/pixels

        shp_polygon_layer = openLayer('vector',shp_polygon, shp_polygon[:-4], "ogr")
        polygon_feature = shp_polygon_layer.getFeatures()
        codetype=code_type(types)
        for featurette in polygon_feature:
            if lsORlsd == "ls":
                results.append([featurette['splitIce'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm,P_average, LAT_ave,LON_ave,Alt_average, pixels,Topo_fc,scal,Input_nuclide,codetype])
            else:
                results.append([featurette['splitIce'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm, LAT_ave,LON_ave,Alt_average, pixels,Topo_fc,scal,Input_nuclide,codetype])
        QgsMapLayerRegistry.instance().removeMapLayer(shp_polygon_layer.id())

    if types == 'GTI':
        if lsORlsd == "ls":
            list_field_names =['Pn_GTI','SFn_GTI', 'Psm_GTI','SFsm_GTI','Pfm_GTI','SFfm_GTI','H_ave', 'LAT_ave','LON_ave','ALT_ave', 'pixels','topo_fc','scaling','nuclide','options']
        else:
            list_field_names =['Pn_GTI','SFn_GTI', 'Psm_GTI','SFsm_GTI','Pfm_GTI','SFfm_GTI','LAT_ave','ALT_ave','LON_ave','pixels','topo_fc','scaling','nuclide','options']
    else:
        if lsORlsd == "ls":
            list_field_names =['Pn_TI','SFn_TI', 'Psm_TI','SFsm_TI','Pfm_TI','SFfm_TI','H_ave', 'LAT_ave','LON_ave','ALT_ave', 'pixels','topo_fc','scaling','nuclide','options']
        else:
            list_field_names =['Pn_TI','SFn_TI', 'Psm_TI','SFsm_TI','Pfm_TI','SFfm_TI','LAT_ave','ALT_ave','LON_ave','pixels','topo_fc','scaling','nuclide','options']
                           
    updateResults(outSHP, list_field_names, results, 'splitIce')

    if use_VDM == True:
        if lsORlsd == "ls":
            list_fields_todelete = ['Psm_P','SFsm_P','Pfm_P','SFfm_P','splitIce']
            if use_geol_mask == True:
                VDM_process(outSHP, VDM, ConcField,SLHLP_Be_neutron, "GTI")
                list_fields_todelete = list_fields_todelete+['splitQz']
            else:
                VDM_process(outSHP, VDM, ConcField,SLHLP_Be_neutron, "TI")
                list_fields_todelete = list_fields_todelete+['split']
        else:
            list_fields_todelete = ['H_ave','splitIce']
            if use_geol_mask == True:
                VDM_process_lsd (outSHP,ConcField,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon, Input_nuclide,"GTI")
                list_fields_todelete = list_fields_todelete+['splitQz']
            else:
                VDM_process_lsd (outSHP,ConcField,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon, Input_nuclide,"TI")
                list_fields_todelete = list_fields_todelete+['split']
    else:
        if lsORlsd == "ls":
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','Time','splitIce']
            if use_geol_mask == True:
                list_fields_todelete = list_fields_todelete+['splitQz']
            else:
                list_fields_todelete = list_fields_todelete+['split']
        else:
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','H_ave','Time','splitIce']
            if use_geol_mask == True:
                list_fields_todelete = list_fields_todelete+['splitQz']
            else:
                list_fields_todelete = list_fields_todelete+['split']
        
    deleteFields(outSHP, list_fields_todelete)

    if use_geol_mask == True:
        deleteFields(inGeolShp, ['splitQz'])
    else:
        deleteFields(shp, ['split'])

    selectSaveByAttribute(outSHP, '\"pixels\">=\'0\'', outGLIMS)
    
    QgsMapLayerRegistry.instance().removeMapLayer(outSHP.id())


def glims_noTS(shp, inDEM, temp_dir, use_geol_mask, inGeolShp,qz_express,outQZshp, use_VDM, VDM, ConcField, inGLIMS, outGLIMS, types,lsORlsd,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide):

    if use_geol_mask == True:    
        temp = temp_dir + '/' + 'geolSHP_temp.shp'
        selectSaveByAttribute(inGeolShp, qz_express, temp)
        temp2 = temp_dir + '/' + 'geolSHP_temp2.shp'
        processing.runalg('qgis:clip', shp, temp, temp2)
        outSHP = openLayer('vector',temp2, temp2[:-4], "ogr")
        selectSaveByAttribute(outSHP,'\"pixels\">=\'0\'',outQZshp)
        inSHP=outQZshp
        inSHP = openLayer('vector', temp2, temp2[:-4], "ogr")
        assignSplitNumber(inSHP, "splitQz")
        list_shp_Qz=splitSHP(inSHP, "splitQz", temp_dir)
    else:
        inSHP=shp
        assignSplitNumber(inSHP, "split")
        list_shp_0=splitSHP(inSHP, "split", temp_dir)

    temp3 = temp_dir + '/' + 'out_ice_temp.shp'
    processing.runalg('qgis:difference',inSHP,inGLIMS, True,temp3)

    outSHP = openLayer('vector', temp3, temp3[:-4], "ogr")

    addFields(outSHP, list_fields(types))

    assignSplitNumber(outSHP, "splitIce")

    results = []

    # Create one file for each polygon of the shapefile   
    list_shp = splitSHP(outSHP, "splitIce", temp_dir)
    length=len(list_shp)
    
    for i in range (0,length):
        shp_polygon=list_shp[i]
        polygonNames = createPolyNames(shp_polygon)

        clipRaster_byLayer(inDEM, shp_polygon, polygonNames[0])

        DEM_loaded = openLayer('raster', polygonNames[0], shp_polygon[:-4] + 'dtm')
        extent = getExtent(DEM_loaded)
        QgsMapLayerRegistry.instance().removeMapLayer(DEM_loaded.id())

        processing.runalg("gdalogr:translate", polygonNames[0], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[2])

        if lsORlsd == "ls": 
            ppg, SUM_Pre, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK = processNumpy(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon)        
        elif lsORlsd == "lsd":
            ppg, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide)
        else:
            ppg, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD_poly(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide)

        if use_geol_mask == True:
            shp_polygon_Qz=list_shp_Qz[i]
            pixels=pixel_count(shp_polygon_Qz,inDEM)       
        else:
            shp_polygon_0=list_shp_0[i]
            pixels=pixel_count(shp_polygon_0,inDEM)

        ProdDEM_Ben = SUM_Pn / pixels
        SF_Ben = SUM_SFn / pixels
        ProdDEM_Besm = SUM_Psm / pixels
        SF_Besm = SUM_SFsm / pixels
        ProdDEM_Befm = SUM_Pfm / pixels
        SF_Befm = SUM_SFfm / pixels
        Alt_average = SUM_DEM / pixels
        if lsORlsd == "ls":
            P_average=SUM_Pre/pixels
            scal=1
        else:
            P_average=0
            scal=0
        LAT_ave=SUM_LAT/pixels
        LON_ave=SUM_LONG/pixels

        shp_polygon_layer = openLayer('vector',shp_polygon, shp_polygon[:-4], "ogr")
        polygon_feature = shp_polygon_layer.getFeatures()

        codetype=code_type(types)
        for featurette in polygon_feature:
            if lsORlsd == "ls":
                results.append([featurette['splitIce'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm,P_average, LAT_ave,LON_ave,Alt_average, pixels,scal,Input_nuclide,codetype])
            else:
                results.append([featurette['splitIce'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm, LAT_ave,LON_ave,Alt_average, pixels,scal,Input_nuclide,codetype])
        QgsMapLayerRegistry.instance().removeMapLayer(shp_polygon_layer.id())

    if types == 'GI':
        if lsORlsd == "ls":
            list_field_names =['Pn_GI','SFn_GI', 'Psm_GI','SFsm_GI','Pfm_GI','SFfm_GI','H_ave', 'LAT_ave','LON_ave','ALT_ave', 'pixels','scaling','nuclide','options']
        else:
            list_field_names =['Pn_GI','SFn_GI', 'Psm_GI','SFsm_GI','Pfm_GI','SFfm_GI','LAT_ave','ALT_ave','LON_ave','pixels','scaling','nuclide','options']
    else:
        if lsORlsd == "ls":
            list_field_names =['Pn_I','SFn_I','Psm_I','SFsm_I','Pfm_I','SFfm_I','H_ave', 'LAT_ave','LON_ave','ALT_ave', 'pixels','scaling','nuclide','options']
        else:
            list_field_names =['Pn_I','SFn_I', 'Psm_I','SFsm_I','Pfm_I','SFfm_I','LAT_ave','ALT_ave','LON_ave', 'pixels','scaling','nuclide','options']
        
    updateResults(outSHP, list_field_names, results, 'splitIce')

    if use_VDM == True:
        if lsORlsd == "ls":      
            list_fields_todelete = ['Psm_P','SFsm_P','Pfm_P','SFfm_P','splitIce']
            if use_geol_mask == True:
                VDM_process(outSHP, VDM, ConcField,SLHLP_Be_neutron, "GI")
                list_fields_todelete = list_fields_todelete+['splitQz']
            else:
                VDM_process(outSHP, VDM, ConcField,SLHLP_Be_neutron, "I")
                list_fields_todelete = list_fields_todelete+['split']
        else:
            list_fields_todelete = ['H_ave','splitIce']
            if use_geol_mask == True:
                VDM_process_lsd (outSHP,ConcField,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon, Input_nuclide,"GI")
                list_fields_todelete = list_fields_todelete+['splitQz']
            else:
                VDM_process_lsd (outSHP,ConcField,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon, Input_nuclide,"I")
                list_fields_todelete = list_fields_todelete+['split']
    else:
        if lsORlsd == "ls":
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','Time','splitIce']
            if use_geol_mask == True:
                list_fields_todelete = list_fields_todelete+['splitQz']
            else:
                list_fields_todelete = list_fields_todelete+['split']
        else:
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','H_ave','Time','splitIce']
            if use_geol_mask == True:
                list_fields_todelete = list_fields_todelete+['splitQz']
            else:
                list_fields_todelete = list_fields_todelete+['split']
        
    deleteFields(outSHP, list_fields_todelete)

    if use_geol_mask == True:
        deleteFields(inGeolShp, ['splitQz'])
    else:
        deleteFields(shp, ['split'])
    
    selectSaveByAttribute(outSHP, '\"pixels\">=\'0\'', outGLIMS)

    QgsMapLayerRegistry.instance().removeMapLayer(outSHP.id())

def pr_Geol_TS(shp,inDEM,outDEM,temp_dir,inGeolShp,qz_express,outQZshp,in_ts,use_VDM,VDM,ConcField,lsORlsd,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide):
    temp = temp_dir + '/' + 'geolSHP_temp.shp'
    selectSaveByAttribute(inGeolShp, qz_express, temp)

    temp2 = temp_dir + '/' + 'geolSHP_temp2.shp'
    processing.runalg('qgis:clip', shp, temp, temp2)

    outSHP = openLayer('vector',temp2, temp2[:-4], "ogr")

    types='GT'
    addFields(outSHP, list_fields('GT'))

    assignSplitNumber(outSHP)

    results = []

    # Create one file for each polygon of the shapefile
    list_shp = splitSHP(outSHP, "split", temp_dir)
    length=len(list_shp)
    
    for i in range (0,length):
        shp_polygon=list_shp[i]
        polygonNames = createPolyNames(shp_polygon)

        clipRaster_byLayer(inDEM, shp_polygon, polygonNames[0])
        clipRaster_byLayer(in_ts, shp_polygon, polygonNames[1])


        DEM_loaded = openLayer('raster',polygonNames[0], shp_polygon[:-4] + 'dtm')
        extent = getExtent(DEM_loaded)
        QgsMapLayerRegistry.instance().removeMapLayer(DEM_loaded.id())

        processing.runalg("gdalogr:translate", polygonNames[0], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[2])
        processing.runalg("gdalogr:translate", polygonNames[1], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[3])

        if lsORlsd == "ls": 
            pp, SUM_Pre, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK = processNumpy(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,polygonNames[3])        
        elif lsORlsd == "lsd":
            pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide,polygonNames[3])
        else:
            pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD_poly(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide,polygonNames[3])


        # Calculation of average production rates in the drainage basins:
        ProdDEM_Ben = SUM_Pn / pp
        SF_Ben = SUM_SFn / pp
        ProdDEM_Besm = SUM_Psm / pp
        SF_Besm = SUM_SFsm / pp
        ProdDEM_Befm = SUM_Pfm / pp
        SF_Befm = SUM_SFfm / pp
        Alt_average = SUM_DEM / pp
        if lsORlsd == "ls":
            P_average=SUM_Pre/pp
            scal=1
        else:
            P_average=0
            scal=0
        LAT_ave=SUM_LAT/pp
        LON_ave=SUM_LONG/pp
        Topo_fc=SUM_MASK/pp

        shp_polygon_layer = openLayer('vector',shp_polygon, shp_polygon[:-4], "ogr")
        polygon_feature = shp_polygon_layer.getFeatures()
        codetype=code_type('GT')
        for featurette in polygon_feature:
            if lsORlsd == "ls":
                results.append([featurette['split'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm,P_average, LAT_ave,LON_ave,Alt_average, pp,Topo_fc,scal,Input_nuclide,codetype])
            else:
                results.append([featurette['split'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm, LAT_ave,LON_ave,Alt_average, pp,Topo_fc,scal,Input_nuclide,codetype])

        QgsMapLayerRegistry.instance().removeMapLayer(shp_polygon_layer.id())

    if lsORlsd == "ls":
        list_field_names =['Pn_GT','SFn_GT', 'Psm_GT','SFsm_GT','Pfm_GT','SFfm_GT','H_ave', 'LAT_ave','LON_ave','ALT_ave', 'pixels','topo_fc','scaling','nuclide','options']
    else:
        list_field_names =['Pn_GT','SFn_GT', 'Psm_GT','SFsm_GT','Pfm_GT','SFfm_GT','LAT_ave','ALT_ave','LON_ave', 'pixels','topo_fc','scaling','nuclide','options']
        
    updateResults(outSHP,list_field_names,results)

    if use_VDM == True:
        if lsORlsd == "ls":
            VDM_process(outSHP, VDM, ConcField,SLHLP_Be_neutron, "GT")
            list_fields_todelete = ['Psm_P','SFsm_P','Pfm_P','SFfm_P','split']
        else:
            VDM_process_lsd (outSHP,ConcField,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon, Input_nuclide,"GT")
            list_fields_todelete = ['H_ave','split']
    else:
        if lsORlsd == "ls":
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','Time','split']
        else:
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','H_ave','Time','split']

    deleteFields(outSHP,list_fields_todelete)

    selectSaveByAttribute(outSHP,'\"pixels\">=\'0\'',outQZshp)

    QgsMapLayerRegistry.instance().removeMapLayer(outSHP.id())

    return pp

def pr_Geol_noTS(shp,inDEM,outDEM,temp_dir,inGeolShp,qz_express,outQZshp,use_VDM,VDM,ConcField,lsORlsd,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide):
    temp = temp_dir + '/' + 'geolSHP_temp.shp'
    selectSaveByAttribute(inGeolShp, qz_express, temp)

    temp2 = temp_dir + '/' + 'geolSHP_temp2.shp'
    processing.runalg('qgis:clip', shp, temp, temp2)

    outSHP = openLayer('vector', temp2, temp2[:-4], "ogr")

    types='G'
    addFields(outSHP, list_fields('G'))

    assignSplitNumber(outSHP)

    results = []

    # Create one file for each polygon of the shapefile
    list_shp = splitSHP(outSHP, "split", temp_dir)

    for shp_polygon in list_shp:

        polygonNames = createPolyNames(shp_polygon)

        clipRaster_byLayer(inDEM, shp_polygon, polygonNames[0])

        DEM_loaded = openLayer('raster', polygonNames[0], shp_polygon[:-4] + 'dtm')
        extent = getExtent(DEM_loaded)
        QgsMapLayerRegistry.instance().removeMapLayer(DEM_loaded.id())

        processing.runalg("gdalogr:translate", polygonNames[0], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[2])

        if lsORlsd == "ls": 
            pp, SUM_Pre, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK = processNumpy(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon)        
        elif lsORlsd == "lsd":
            pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide)
        else:
            pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD_poly(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide)

        # Calculation of average production rates in the drainage basins:
        ProdDEM_Ben = SUM_Pn / pp
        SF_Ben = SUM_SFn / pp
        ProdDEM_Besm = SUM_Psm / pp
        SF_Besm = SUM_SFsm / pp
        ProdDEM_Befm = SUM_Pfm / pp
        SF_Befm = SUM_SFfm / pp
        Alt_average = SUM_DEM / pp
        if lsORlsd == "ls":
            P_average=SUM_Pre/pp
            scal=1
        else:
            P_average=0
            scal=0
        LAT_ave=SUM_LAT/pp
        LON_ave=SUM_LONG/pp

        shp_polygon_layer = openLayer('vector',shp_polygon, shp_polygon[:-4], "ogr")
        polygon_feature = shp_polygon_layer.getFeatures()
        codetype=code_type('G')
        for featurette in polygon_feature:
            if lsORlsd == "ls":
                results.append([featurette['split'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm,P_average, LAT_ave,LON_ave,Alt_average, pp,scal,Input_nuclide,codetype])
            else:
                results.append([featurette['split'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm, LAT_ave,LON_ave,Alt_average, pp,scal,Input_nuclide,codetype])
        QgsMapLayerRegistry.instance().removeMapLayer(shp_polygon_layer.id())

    if lsORlsd == "ls":
        list_field_names =['Pn_G','SFn_G', 'Psm_G','SFsm_G','Pfm_G','SFfm_G','H_ave', 'LAT_ave','LON_ave','ALT_ave', 'pixels','scaling','nuclide','options']
    else:
        list_field_names =['Pn_G','SFn_G', 'Psm_G','SFsm_G','Pfm_G','SFfm_G','LAT_ave','LON_ave','ALT_ave', 'pixels','scaling','nuclide','options']
  
    updateResults(outSHP, list_field_names, results)

    if use_VDM == True:
        if lsORlsd == "ls":
            VDM_process(outSHP, VDM, ConcField,SLHLP_Be_neutron, "G")
            list_fields_todelete = ['Psm_P','SFsm_P','Pfm_P','SFfm_P','split']
        else:
            VDM_process_lsd (outSHP,ConcField,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon, Input_nuclide,"G")
            list_fields_todelete = ['H_ave','split']
    else:
        if lsORlsd == "ls":
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','Time','split']
        else:
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','H_ave','Time','split']

    deleteFields(outSHP, list_fields_todelete)

    selectSaveByAttribute(outSHP, '\"pixels\">=\'0\'', outQZshp)

    QgsMapLayerRegistry.instance().removeMapLayer(outSHP.id())

    return pp

def pr_noGeol_TS(shp,inDEM,outDEM,temp_dir,in_ts,use_VDM,VDM,ConcField,lsORlsd,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide):

    addFields(shp, list_fields('T'))

    assignSplitNumber(shp)

    results = []

    list_shp = splitSHP(shp, "split", temp_dir)

    for shp_polygon in list_shp:

        polygonNames = createPolyNames(shp_polygon)

        clipRaster_byLayer(inDEM, shp_polygon, polygonNames[0])
        clipRaster_byLayer(in_ts, shp_polygon, polygonNames[1])

        DEM_loaded = openLayer('raster',polygonNames[0], shp_polygon[:-4] + 'dtm')
        extent = getExtent(DEM_loaded)
        QgsMapLayerRegistry.instance().removeMapLayer(DEM_loaded.id())

        processing.runalg("gdalogr:translate", polygonNames[0], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[2])
        processing.runalg("gdalogr:translate", polygonNames[1], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[3])

        if lsORlsd == "ls":
            pp, SUM_Pre, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK = processNumpy(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,polygonNames[3])        
        elif lsORlsd == "lsd":
            pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide,polygonNames[3])
        else:
            pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD_poly(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide,polygonNames[3])

        
        # Calculation of average production rates in the drainage basins:
        ProdDEM_Ben = SUM_Pn / pp
        SF_Ben = SUM_SFn / pp
        ProdDEM_Besm = SUM_Psm / pp
        SF_Besm = SUM_SFsm / pp
        ProdDEM_Befm = SUM_Pfm / pp
        SF_Befm = SUM_SFfm / pp
        Alt_average = SUM_DEM / pp
        if lsORlsd == "ls":
            P_average=SUM_Pre/pp
            scal=1
        else:
            P_average=0
            scal=0
        LAT_ave=SUM_LAT/pp
        LON_ave=SUM_LONG/pp
        Topo_fc=SUM_MASK/pp

        shp_polygon_layer = openLayer('vector',shp_polygon, shp_polygon[:-4], "ogr")

        polygon_feature = shp_polygon_layer.getFeatures()
        codetype=code_type('T')
        for featurette in polygon_feature:
            if lsORlsd == "ls":
                results.append([featurette['split'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm,P_average, LAT_ave,LON_ave,Alt_average, pp,Topo_fc,scal,Input_nuclide,codetype])
            else:
                results.append([featurette['split'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm, LAT_ave,LON_ave,Alt_average, pp,Topo_fc,scal,Input_nuclide,codetype])
                
        QgsMapLayerRegistry.instance().removeMapLayer(shp_polygon_layer.id())

    if lsORlsd == "ls":
        list_field_names =['Pn_T','SFn_T', 'Psm_T','SFsm_T','Pfm_T','SFfm_T','H_ave', 'LAT_ave','LON_ave','ALT_ave', 'pixels','topo_fc','scaling','nuclide','options']
    else:
        list_field_names =['Pn_T','SFn_T', 'Psm_T','SFsm_T','Pfm_T','SFfm_T','LAT_ave','LON_ave','ALT_ave', 'pixels','topo_fc','scaling','nuclide','options']
        
    updateResults(shp, list_field_names, results)

    if use_VDM == True:
        if lsORlsd == "ls":
            VDM_process(shp, VDM, ConcField,SLHLP_Be_neutron, "T")
            list_fields_todelete = ['Psm_P','SFsm_P','Pfm_P','SFfm_P','split']
        else:
            VDM_process_lsd (shp,ConcField,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon, Input_nuclide,"T")
            list_fields_todelete = ['H_ave','split']
    else:
        if lsORlsd == "ls":
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','Time','split']
        else:
            list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','H_ave','Time','split']    

    deleteFields(shp, list_fields_todelete)
    
    return pp
#
# def pr_noGeol_noTS(shp,inDEM,outDEM,temp_dir,use_VDM,VDM,ConcField,lsORlsd,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide):
#
#     addFields(shp, list_fields('None'))
#
#     assignSplitNumber(shp)
#
#     results = []
#
#     # Create one file for each polygon of the shapefile
#     list_shp = splitSHP(shp, "split", temp_dir)
#
#     for shp_polygon in list_shp:
#
#         polygonNames = createPolyNames(shp_polygon)
#
#         clipRaster_byLayer(inDEM, shp_polygon, polygonNames[0])
#
#         DEM_loaded = openLayer('raster', polygonNames[0], shp_polygon[:-4] + 'dtm')
#         extent = getExtent(DEM_loaded)
#         QgsMapLayerRegistry.instance().removeMapLayer(DEM_loaded.id())
#
#         processing.runalg("gdalogr:translate", polygonNames[0], 100, True, "-9999", 0, "", extent, False, 5, 4, 75, 6, 1, False, 0, False, "", polygonNames[2])
#
#         if lsORlsd == "ls":
#             pp, SUM_Pre, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon)
#         elif lsORlsd == "lsd":
#             pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide)
#         else:
#             pp, SUM_LAT, SUM_LONG, SUM_DEM, SUM_Pn, SUM_Psm, SUM_Pfm, SUM_SFn, SUM_SFsm, SUM_SFfm, SUM_MASK= processNumpy_LSD_poly(polygonNames[2],SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,Input_nuclide)
#
#
#         # Calculation of average production rates in the drainage basins:
#         ProdDEM_Ben = SUM_Pn / pp
#         SF_Ben = SUM_SFn / pp
#         ProdDEM_Besm = SUM_Psm / pp
#         SF_Besm = SUM_SFsm / pp
#         ProdDEM_Befm = SUM_Pfm / pp
#         SF_Befm = SUM_SFfm / pp
#         Alt_average = SUM_DEM / pp
#         LAT_ave=SUM_LAT/pp
#         LON_ave=SUM_LONG/pp
#         if lsORlsd == "ls":
#             P_average=SUM_Pre/pp
#             scal=1
#         else:
#             P_average=0
#             scal=0
#         shp_polygon_layer = openLayer('vector', shp_polygon, shp_polygon[:-4], "ogr")
#
#         polygon_feature = shp_polygon_layer.getFeatures()
#         codetype=code_type('None')
#         for featurette in polygon_feature:
#             if lsORlsd == "ls":
#                 results.append([featurette['split'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm,P_average, LAT_ave,LON_ave,Alt_average, pp,scal,Input_nuclide,codetype])
#             else:
#                 results.append([featurette['split'],ProdDEM_Ben, SF_Ben, ProdDEM_Besm, SF_Besm, ProdDEM_Befm, SF_Befm, LAT_ave,LON_ave,Alt_average, pp,scal,Input_nuclide,codetype])
#
#         QgsMapLayerRegistry.instance().removeMapLayer(shp_polygon_layer.id())
#
#     if lsORlsd == "ls":
#         list_field_names =['Pn_','SFn_', 'Psm_','SFsm_','Pfm_','SFfm_','H_ave', 'LAT_ave','LON_ave','ALT_ave', 'pixels','scaling','nuclide','options']
#     else:
#         list_field_names =['Pn_','SFn_', 'Psm_','SFsm_','Pfm_','SFfm_','LAT_ave','LON_ave','ALT_ave', 'pixels','scaling','nuclide','options']
#
#     updateResults(shp, list_field_names, results)
#
#     if use_VDM == True:
#         if lsORlsd == "ls":
#             VDM_process(shp, VDM, ConcField,SLHLP_Be_neutron, "None")
#             list_fields_todelete = ['Psm_P','SFsm_P','Pfm_P','SFfm_P','split']
#         else:
#             VDM_process_lsd (shp,ConcField,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon, Input_nuclide,"None")
#             list_fields_todelete = ['H_ave','split']
#     else:
#         if lsORlsd == "ls":
#             list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','Time','split']
#         else:
#             list_fields_todelete = ['Pn_P','SFn_P','Psm_P','SFsm_P','Pfm_P','SFfm_P','H_ave','Time','split']
#
#     deleteFields(shp, list_fields_todelete)
#
#     return pp
#
def pr_noGeol_noTS(shp, inDEM, outDEM, temp_dir, use_VDM, VDM, ConcField, Input_scaling_scheme,
                   SLHLP_Be_neutron, SLHLP_Be_fastmuon, SLHLP_Be_slowmuon, nuclide):
    """
    Process the data without geological or timescale mask.

    Args:
        shp: Input shapefile.
        inDEM: Input DEM.
        outDEM: Output DEM.
        temp_dir: Temporary directory for processing.
        use_VDM: Use vertical density mask.
        VDM: Vertical density mask data.
        ConcField: Concentration field.
        Input_scaling_scheme: Scaling scheme input.
        SLHLP_Be_neutron, SLHLP_Be_fastmuon, SLHLP_Be_slowmuon: Constants for scaling.
        nuclide: Nuclide being processed.

    Returns:
        Processed result (pp).
    """
    pp = None  # Default value for pp

    try:
        # Core processing logic here
        if some_condition:  # Replace with actual condition
            pp = some_processing_function(shp, inDEM, outDEM)
        else:
            print("Error: Missing or invalid inputs for processing.")
            return None

    except Exception as e:
        print(f"Error in pr_noGeol_noTS: {e}")
        return None

    return pp

#===========================================================================================================================================================================================
#Main core of the Basinga program: define the SLHL and call the main functions accordingly to the option selected.
#===========================================================================================================================================================================================
def pr_processes(shp,inDEM,outDEM,temp_dir,Input_nuclide, use_geol_mask,inGeolShp,qz_express,outQZshp,use_ts,in_ts,use_VDM,ConcField,use_glims_mask,inGLIMS,outGLIMS,Input_scaling_scheme):
    # main function, lots of parameters, direct to the right process depending on the desired corrections

    # Output the DEM of the basin
    temp = temp_dir + '/' + 'mergedSHP_temp.shp'
    # processing.runalg('qgis:dissolve', shp, True, None, temp)
    ##MZ
    processing.run("native:dissolve", {
        'INPUT': shp,
        'FIELD': [],  # Leave empty for dissolving all features
        'OUTPUT': temp
    })
    clipRaster_byLayer(inDEM, temp, outDEM)

    if Input_nuclide == "10Be":
        if Input_scaling_scheme == "ls":
            ##Lal/Stone + ERA40 atmosphere + Muscheler:
            SLHLP_Be=4.18 #Table 7 Martin et al. (2017) : cREP
            nuclide=10
        else:
            ##LSD + ERA40 atmosphere + Muscheler:
            SLHLP_Be=4.14 #Table 7 Martin et al. (2017) : cREP
            nuclide=10

    elif Input_nuclide == "3He":
    ##Please note that the code was originally written for 10Be only, for simplicty we keep the same for the parameters with "_Be" but it does consider different nuclide
        if Input_scaling_scheme == "ls":
            ##Lal/Stone + ERA40 atmosphere + Muscheler:
            SLHLP_Be=122 #Table 7Martin et al. (2017) : cREP 
            nuclide=3
        else: 
            ##LSD + ERA40 atmosphere + Muscheler:
            SLHLP_Be=125 #Table 7 Martin et al. (2017) : cREP
            nuclide=3
        
    elif Input_nuclide == "21Ne":
    ##Please note that the code was originally written for 10Be only, for simplicty we keep the same for the parameters with "_Be" but it does consider different nuclide
        if Input_scaling_scheme == "ls":
            ##Lal/Stone + ERA40 atmosphere + Muscheler:
            SLHLP_Be=4.18/4.12 #10Be/21Ne=4.12 +/- 0.17 see Kober et al., EPSL 2011 and Table 7 Martin et al. (2017) : cREP
            nuclide=21
        else: 
            #LSD + ERA40 atmosphere + Muscheler:
            SLHLP_Be=4.14/4.12 #10Be/21Ne=4.12 +/- 0.17 see Kober et al., EPSL 2011 and Table 7 Martin et al. (2017) : cREP
            nuclide=21
        
    SLHLP_Be_neutron = fsp*SLHLP_Be
    SLHLP_Be_fastmuon=fmu_fastmuon*SLHLP_Be
    SLHLP_Be_slowmuon=fmu_slowmuon*SLHLP_Be

    # direct to the right process depending on the desired corrections
    if use_geol_mask == True:
        if use_ts == True:
            if use_glims_mask == False:
                pp=pr_Geol_TS(shp, inDEM, outDEM, temp_dir, inGeolShp, qz_express, outQZshp, in_ts,use_VDM,VDM,ConcField,Input_scaling_scheme,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide)
            elif use_glims_mask == True:
                glims_TS(shp,inDEM,temp_dir,in_ts,use_geol_mask, inGeolShp,qz_express,outQZshp, use_VDM,VDM,ConcField,inGLIMS,outGLIMS,'GTI',Input_scaling_scheme,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide)
        else:
            if use_glims_mask == False:
                pp=pr_Geol_noTS(shp, inDEM, outDEM, temp_dir, inGeolShp, qz_express, outQZshp,use_VDM,VDM,ConcField,Input_scaling_scheme,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide)
            elif use_glims_mask == True:
                glims_noTS(shp,inDEM,temp_dir,use_geol_mask, inGeolShp,qz_express,outQZshp, use_VDM,VDM,ConcField,inGLIMS,outGLIMS,'GI',Input_scaling_scheme,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide)
    else:
        if use_ts == True:
            if use_glims_mask == False:
                pp=pr_noGeol_TS(shp, inDEM, outDEM, temp_dir, in_ts,use_VDM,VDM,ConcField,Input_scaling_scheme,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide)
            elif use_glims_mask == True:
                glims_TS(shp,inDEM,temp_dir,in_ts,use_geol_mask, inGeolShp,qz_express,outQZshp, use_VDM,VDM,ConcField,inGLIMS,outGLIMS,'TI',Input_scaling_scheme,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide)
        else :
            if use_glims_mask == False:
                pp=pr_noGeol_noTS(shp, inDEM, outDEM, temp_dir ,use_VDM,VDM,ConcField,Input_scaling_scheme,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide)
            elif use_glims_mask == True:
                glims_noTS(shp,inDEM,temp_dir,use_geol_mask, inGeolShp,qz_express,outQZshp, use_VDM,VDM,ConcField,inGLIMS,outGLIMS,'I',Input_scaling_scheme,SLHLP_Be_neutron,SLHLP_Be_fastmuon,SLHLP_Be_slowmuon,nuclide)

