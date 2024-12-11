# -*- coding: utf-8 -*-
# ==================================================================
#===========================================================================================================================================================================================
#Auxiliary functions needed for the calculations
#===========================================================================================================================================================================================

# Basic mathematical functions
import numpy as np
import math
from math import exp
from scipy import interpolate

# Import parameters for calculation
from parameters_definition import *


def list_fields(type):
    # Return a list of the names of the fields used in the process, requires a string specifying the type of corrections used : "GTI", "TI", "T"...
    # G stands for geological correction
    # T for Topographic shielding correction
    # I for Ice cover

    # Each field name is followed by "Double" which specify the kind of value that can be stored in the field.
    if type == 'GTI':
        list_new_fields = [["Pn_GTI","Double"],["SFn_GTI","Double"], ["Psm_GTI","Double"], ["SFsm_GTI","Double"], ["Pfm_GTI","Double"], ["SFfm_GTI","Double"],["topo_fc","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["LON_ave","Double"],["ALT_ave","Double"], ["pixels","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type == 'GI':
        list_new_fields = [["Pn_GI","Double"],["SFn_GI","Double"], ["Psm_GI","Double"], ["SFsm_GI","Double"], ["Pfm_GI","Double"], ["SFfm_GI","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"], ["LON_ave","Double"],["pixels","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type == 'TI':
        list_new_fields = [["Pn_TI","Double"],["SFn_TI","Double"], ["Psm_TI","Double"], ["SFsm_TI","Double"], ["Pfm_TI","Double"], ["SFfm_TI","Double"],["topo_fc","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"],["LON_ave","Double"], ["pixels","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type == 'I':
        list_new_fields = [["Pn_I","Double"],["SFn_I","Double"], ["Psm_I","Double"], ["SFsm_I","Double"], ["Pfm_I","Double"], ["SFfm_I","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"], ["LON_ave","Double"],["pixels","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type == 'GT':
        list_new_fields = [["Pn_GT","Double"],["SFn_GT","Double"], ["Psm_GT","Double"], ["SFsm_GT","Double"], ["Pfm_GT","Double"], ["SFfm_GT","Double"],["topo_fc","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"],["LON_ave","Double"], ["pixels","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type == 'G':
        list_new_fields = [["Pn_G","Double"],["SFn_G","Double"], ["Psm_G","Double"], ["SFsm_G","Double"], ["Pfm_G","Double"], ["SFfm_G","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"], ["LON_ave","Double"],["pixels","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    elif type == 'T':
        list_new_fields = [["Pn_T","Double"],["SFn_T","Double"], ["Psm_T","Double"], ["SFsm_T","Double"], ["Pfm_T","Double"], ["SFfm_T","Double"],["topo_fc","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"],["LON_ave","Double"], ["pixels","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]

    else :
        list_new_fields = [["Pn_","Double"],["SFn_","Double"], ["Psm_","Double"], ["SFsm_","Double"], ["Pfm_","Double"], ["SFfm_","Double"]]
        list_new_fields = list_new_fields + [["H_ave","Double"], ["LAT_ave","Double"], ["ALT_ave","Double"],["LON_ave","Double"], ["pixels","Double"]]
        list_new_fields = list_new_fields + [["Pn_P","Double"],["SFn_P","Double"], ["Psm_P","Double"], ["SFsm_P","Double"], ["Pfm_P","Double"], ["SFfm_P","Double"],["Time", "Double"]]
    
    return list_new_fields

def checkFields(inSHP,listFields):
    #check wether or not the inSHP file already includes the fields needed to store the results of the requested calculations  as function of the selected options
    list_field_u=[field.name() for field in inSHP.pendingFields()]
    list_field0=[x.encode('UTF8') for x in list_field_u]
    result=0
    for field in list_field0:
        ref=listFields[0]
        ref=ref[0]
        if field == ref:
            result=1
    return result
    
def addFields(inSHP,listFields):
    # Add new field to a vector layer
    # inSHP must be an open layer, listFields a list of the new field names and their type : Int or Double
    # listFields = [[name,type],["newfield","Int"]]
    caps = inSHP.dataProvider().capabilities()
    if caps & QgsVectorDataProvider.AddFeatures:
        for field in listFields: # add each field specified, the type is respected
            if field[1] == "Int":
                inSHP.dataProvider().addAttributes([QgsField(field[0], QVariant.Int)])
            elif field[1] == "Double":
                inSHP.dataProvider().addAttributes([QgsField(field[0], QVariant.Double)])

def assignSplitNumber (inSHP,field="split",type = "Int"):
    # assign a split number by creating a split flied and counting the polygons/features, inSHP must be an open layer, a field name and a value type can be specified

    addFields(inSHP,[[field, type]])

    iter = inSHP.getFeatures() # get a list of the feature/polygon objects
    i = 1
    inSHP.startEditing() # open the editor mode
    idx = inSHP.fieldNameIndex(field) # get the index of the field "field"
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

def openLayer(type,path,legendName,pilot='ogr'):
    #open vector or raster layer, requires the type ("raster" or "layer"), the path of the wanted file, the name that will be used
    # in the qgis legend and the pilot for vector file (defautl set as "ogr", no change needed)
    # Return the opened layer as an qgis layer object
    if type == 'raster':
        objectName = QgsRasterLayer(path, legendName)
    elif type == 'vector':
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
    processing.runalg('qgis:splitvectorlayer', inSHP, fieldName, dir) # qgis command to split a vector layer to x files containing each one the polygons sharing a common field value

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
    v = qgis.utils.QGis.QGIS_VERSION
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
        idx = inSHP.fieldNameIndex(field)
        list_fields_index.append(idx)

    return list_fields_index

def updateResults(inSHP,listFieldName,results,split='split'):
    # requires an opened layer, the list of the field names used, the matrix of the results, and the name of the field used to count the polygons
    list_fields_index = getFieldIndex(inSHP, listFieldName) # get the list the list of the field indexes to work with it
    
    iter = inSHP.getFeatures()
    inSHP.startEditing()

    # copy the values in results in the matching place in the attribute table of inSHP
    i = 0
    for features in iter:
        num = features[split]
        t = -1
        it = 0
        while t == -1:
            num2 = results[it][0]
            if num == num2: # Make sure that the split number match to that the two table present the same line when copying
                t = it
            it += 1
        for j in range(len(list_fields_index)):
            inSHP.changeAttributeValue(features.id(), list_fields_index[j], float(results[t][j + 1])) # we save the results in the wanted fields by specifying the feature, the index of the field and the value
        i += 1

    inSHP.commitChanges()

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

    Input: latitude array and cellsize

    Output: weight

    """
    Rt=6370000 #Earth radius in m
    non_zero = latitude > 0
    lat_min=np.min(latitude[non_zero])
    lat_rad=latitude[non_zero]*pi/180
    r_min=Rt*np.cos(lat_min)
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
