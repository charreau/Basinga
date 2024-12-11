
# -*- coding: utf-8 -*-
#=================================================================================================
#Script which calculates the average erosion

# Basic mathematical functions
import math
from math import sqrt

# Import cosmogenic parameters and the addfield function
from parameters_definition import *
from  PR_processes import addFields

#Uncertainties on the 10Be cosmogenic production rates:
#U_spal=0.09 #Balco et al. (2008) - we used instead the uncertianties given in Table 7 of Martin et al. (2017) cf function "nuclideslected" below
U_sm=0.2
U_fm=0.2

def checkFields_denu(inSHP,ifVDM):
    #check wether or not the inSHP file already includes the fields needed to store the results of the requested calculations  as function of the selected options
    list_field_u=[field.name() for field in inSHP.pendingFields()]
    list_field0=[x.encode('UTF8') for x in list_field_u]
    result=0
    if ifVDM == True:
        list_check=[["denuda_P","Double"],["Uncero_P","Double"]]
        for field in list_field0:
            ref=list_check[0]
            ref=ref[0]
            if field == ref:
                result=1
    else:
        list_check=[["denuda_","Double"],["Uncero_","Double"]]
        for field in list_field0:
            ref=list_check[0]
            ref=ref[0]
            if field == ref:
                result=1
                
    return result

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
        list_fields = ['Pn_GTI','Pfm_GTI','Psm_GTI']        
    elif code ==2:
        list_fields = ['Pn_GI','Pfm_GI','Psm_GI']
    elif code ==3:
        list_fields = ['Pn_TI','Pfm_TI','Psm_TI']
    elif code ==4:
        list_fields = ['Pn_I','Pfm_I','Psm_I']
    elif code ==5:
        list_fields = ['Pn_GT','Pfm_GT','Psm_GT']
    elif code == 6:
        list_fields = ['Pn_G','Pfm_G','Psm_G']
    elif code == 7:
        list_fields = ['Pn_T','Pfm_T','Psm_T']
    else:
        list_fields = ['Pn_','Pfm_','Psm_']

    return list_fields

def optionsSelected_VDM(code):
    if code == 1:
        list_fields = ['Pn_P','Pfm_GTI','Psm_GTI']        
    elif code ==2:
        list_fields = ['Pn_P','Pfm_GI','Psm_GI']
    elif code ==3:
        list_fields = ['Pn_P','Pfm_TI','Psm_TI']
    elif code ==4:
        list_fields = ['Pn_P','Pfm_I','Psm_I']
    elif code ==5:
        list_fields = ['Pn_P','Pfm_GT','Psm_GT']
    elif code == 6:
        list_fields = ['Pn_P','Pfm_G','Psm_G']
    elif code == 7:
        list_fields = ['Pn_P','Pfm_T','Psm_T']
    else:
        list_fields = ['Pn_P','Pfm_','Psm_']

    return list_fields

def denudation(inSHP,density,ConcField,UncConc,ifVDM):
    # Process and return the erosion and its uncertainties of a drainage basin. Requires an opened layer, a average rock density,
    # the field of the found isotopic concentration, the field specifying the uncertainties and a boolean which stipulates the
    # use of paleomagnetic correction during the production rate process

    # Two new fields are created, they will eventually contain the results
    result=checkFields_denu(inSHP,ifVDM)

    if result==0:
        if ifVDM == True:
            liste_new_fields = [["denuda_P","Double"],["Uncero_P","Double"]]
        else:
            liste_new_fields = [["denuda_", "Double"], ["Uncero_", "Double"]]   

        addFields(inSHP,liste_new_fields) # Add the new fields to the pre-existing ones, only if they do not already exist

    # We update the layer
    inSHP.startEditing()
    inSHP.commitChanges()

    #retrieving the setting, i.e. scaling, nuclide and options selected for the production calculation
    list_field_setting=['scaling','nuclide','options']
    list_index_setting = []
    for field in list_field_setting:
        list_index_setting.append(inSHP.fieldNameIndex(field))

    features = inSHP.getFeatures() # Get a lit of the features
    inSHP.startEditing() # Switch to editor mode
    for feature in features:
        scaling=feature.attributes()[list_index_setting[0]]
        nuclide=feature.attributes()[list_index_setting[1]]
        code=feature.attributes()[list_index_setting[2]]
    
    # We specify a list of field used in the process
    if ifVDM == True:
        list_fields_in = optionsSelected_VDM(code)+[ConcField, UncConc,"denuda_P", "Uncero_P"]
    else:
        list_fields_in = optionsSelected(code)+[ConcField, UncConc, "denuda_", "Uncero_"]

    # We get the indexes of the specified fields because we will need them to get and change values from those fields
    list_index_in = []
    for field in list_fields_in:
        list_index_in.append(inSHP.fieldNameIndex(field))

    # We start working on every polygon/feature in the layer
    features = inSHP.getFeatures() # Get a lit of the features
    inSHP.startEditing() # Switch to editor mode
    for feature in features:
        # We get the values of the wanted fields
        # The values are accessed by specifying a feature and the index of the field value = feature.attributes()[index of field]
        ProdR_spal = feature.attributes()[list_index_in[0]]
        ProdR_sm = feature.attributes()[list_index_in[1]]
        ProdR_fm = feature.attributes()[list_index_in[2]]
        Conc = feature.attributes()[list_index_in[3]]
        Conc=float(Conc)
        U_Conc = feature.attributes()[list_index_in[4]]
        U_Conc=float(U_Conc)

        # We start calculation
        U_spal=nuclideSelected(nuclide,scaling)
        UncProdR_spal = ProdR_spal * U_spal
        UncProdR_sm = ProdR_sm * U_sm
        UncProdR_fm = ProdR_fm * U_fm

        P_123 = (ProdR_spal*lght_att_n + ProdR_sm*lght_att_sm + ProdR_fm*lght_att_fm) / density
        Erosion = P_123 / Conc  # cm/yr
        
        dc_c2 = U_Conc / (Conc ** 2)
        A = (dc_c2 * P_123) ** 2
        B = (UncProdR_spal * lght_att_n / (Conc * density)) ** 2
        C = (UncProdR_sm * lght_att_sm / (Conc * density)) ** 2
        D = (UncProdR_fm * lght_att_fm / (Conc * density)) ** 2
        Uncert_ero = sqrt(A + B + C + D)

        # we save the results in the wanted fields by specifying the feature, the index of the field and the value
        inSHP.changeAttributeValue(feature.id(), list_index_in[5], Erosion)
        inSHP.changeAttributeValue(feature.id(), list_index_in[6], Uncert_ero)

    inSHP.commitChanges() # Save the changes

