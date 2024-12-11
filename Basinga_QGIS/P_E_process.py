# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Basinga
                                 A QGIS plugin
 Caculate 10Be production rates and erosion
                              -------------------
        begin                : 2016-06-21
        git sha              : $Format:%H$
        copyright            : (C) 2016 by CRPG - Nancy France
        email                : charreau@crpg.cnrs-nancy.fr
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt5.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from PyQt5.QtGui import QIcon
from PyQt5.QtWidgets import QAction
from qgis.gui import QgsMessageBar
# Initialize Qt resources from file resources.py
from . import resources
# Import the code for the dialog
from .P_E_process_dialog import BasingaDialog, SQLDialog, HelpDialog
from .PR_processes import pr_processes
from .Denudation_process import denudation
import os
import os.path
import stat
import datetime
import sys

from qgis.core import Qgis
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QFileDialog


class Basinga:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(os.getcwd())
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'Basinga_{}.qm'.format(locale))

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        # Create the dialog (after translation) and keep reference
        print("************************  ############### opening BasinDialog ***********************  ###############")
        self.dlg = BasingaDialog()
        self.dlgsql = SQLDialog()
        self.dlghelp = HelpDialog()

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&Basinga')

        self.toolbar = self.iface.addToolBar(u'Basinga')
        self.toolbar.setObjectName(u'Basinga')


        # Define the buttons of the GUI and specify their actions
        #  production rates tab

        # For OK, CANCEL and HELP buttons
        self.dlg.pushButton_cancel_P.clicked.connect(self.cancel)
        self.dlg.pushButton_ok_P.clicked.connect(self.ok_P)
        self.dlg.pushButton_help_P.clicked.connect(self.help)

        # For output field box, first we clear the text box then we define the save method
        self.dlg.lineEdit_output_DEM_basin.clear()
        self.dlg.pushButton_output_DEM_basin.clicked.connect(self.select_output_DEM_basin_file)

        # For output field box, first we clear the text box then we define the save method
        self.dlg.lineEdit_output_qz_polygons.clear()
        self.dlg.pushButton_output_qz_polygons.clicked.connect(self.select_output_qz_polygons)

        # Commented lines, used to grey-out part of the GUI
        # self.dlg.lineEdit_output_qz_polygons.setEnabled(False)
        # self.dlg.pushButton_output_qz_polygons.setEnabled(False)

        # For output field box, first we clear the text box then we define the save method
        self.dlg.lineEdit_output_glims_polygons.clear()
        self.dlg.pushButton_output_glims_polygons.clicked.connect(self.select_output_glims_polygons)

        # Commented lines, used to grey-out part of the GUI
        # self.dlg.lineEdit_output_glims_polygons.setEnabled(False)
        # self.dlg.pushButton_output_glims_polygons.setEnabled(False)

        # For output expression box, first we clear the text box then we define the save method
        self.dlg.lineEdit_express_qz_extract.clear()
        self.dlg.pushButton_sql.clicked.connect(self.sql_expression)

        # Define the action used for the checkboxes
        self.dlg.checkBox_use_geol_mask.clicked.connect(self.check_geol)
        self.dlg.checkBox_use_ts_mask.clicked.connect(self.check_ts)
        self.dlg.checkBox_VDM.clicked.connect(self.check_VDM)
        self.dlg.checkBox_use_glims_mask.clicked.connect(self.check_glims)

        # set booleans for the use of the different corrections
        self.use_geol_mask = False
        self.use_ts_mask = False
        self.use_VDM = False
        self.use_glims_mask = False

        # Commented lines, used to grey-out part of the GUI 
        #if checkBox_use_geol_mask='True'
        #    self.dlg.lineEdit_express_qz_extract.setEnabled(True)
        #    self.dlg.pushButton_sql.setEnabled(True)

        # define radiobutton for ls or lsd algorithm
        self.dlg.radioButton_ls.setChecked(True)
        self.dlg.radioButton_lsd.setChecked(False)
        self.dlg.radioButton_lsdp.setChecked(False)
        self.dlg.radioButton_ls.clicked.connect(self.pickAlgoLS)
        self.dlg.radioButton_lsd.clicked.connect(self.pickAlgoLSD)
        self.dlg.radioButton_lsdp.clicked.connect(self.pickAlgoLSDp)
        self.lsORlsd = 'ls' # set the choice  for algorithm in string format


        #Denudation tab

        # For OK, CANCEL and HELP buttons
        self.dlg.pushButton_cancel_D.clicked.connect(self.cancel)
        self.dlg.pushButton_ok_D.clicked.connect(self.ok_D)
        self.dlg.pushButton_help_D.clicked.connect(self.help)

        # For density box, first we clear the text box then we define the save method
        self.dlg.lineEdit_density.clear()
        self.dlg.checkBox_ifVDM.clicked.connect(self.check_ifVDM_used)
        self.ifVDM = False

        # For expression window

        # For output expression box, first we clear the text box then we define the save method
        self.dlgsql.pushButton_cancel.clicked.connect(self.cancel_sql)
        self.dlgsql.pushButton_ok.clicked.connect(self.ok_sql)

        self.dlgsql.lineEdit_value.clear()


    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('Basinga', message)


    def add_action(
        self,
        icon_path,
        text,
        callback,
        enabled_flag=True,
        add_to_menu=True,
        add_to_toolbar=True,
        status_tip=None,
        whats_this=None,
        parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/Basinga/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'Basinga'),
            callback=self.run,
            parent=self.iface.mainWindow())

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&Basinga'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the toolbar
        del self.toolbar

    def run(self):
        """Run method that performs all the real work"""
        from qgis.core import QgsProject
        # get the list of the opened layers and then their names
        #self.layers = QgsProject.instance().mapLayers().values()
        self.layers = list(QgsProject.instance().mapLayers().values())
        self.layer_list = []
        for layer in self.layers:
            self.layer_list.append(layer.name())

        # if there's no opened layer, output a warning message
        if len(self.layers) == 0:
            self.iface.messageBar().pushMessage("Hello", "No loaded layer, please open at least the shapefile of your basin", level=Qgis.Warning, duration=10)
            return

        # clear the comboBoxes
        self.dlg.comboBox_shp_basin.clear()
        self.dlg.comboBox_DEM_raster.clear()
        self.dlg.comboBox_shp_geol.clear()
        self.dlg.comboBox_ts_raster.clear()
        self.dlg.comboBox_shp_glims.clear()

        # add the list of layers to the comboboxes lists and perform an update if the selected basin shp is changed
        self.dlg.comboBox_shp_basin.addItems(self.layer_list)
        self.dlg.comboBox_shp_basin.currentIndexChanged[str].connect(self.shp_changed)

        # add the list of layers to the comboboxes lists
        self.dlg.comboBox_DEM_raster.addItems(self.layer_list)
        self.dlg.comboBox_shp_geol.addItems(self.layer_list)
        self.dlg.comboBox_ts_raster.addItems(self.layer_list)
        self.dlg.comboBox_shp_glims.addItems(self.layer_list)

        # create the vdm list and add it to the combobox
        #self.VDM_list = ["GLOPIS","MUSCHELER","VDM2000"]
        #self.dlg.comboBox_VDM.clear()
        #self.dlg.comboBox_VDM.addItems(self.VDM_list)

        # create the nuclide list and add it to the combobox
        self.nuclide_list = ["10Be","3He","21Ne"]
        self.dlg.comboBox_nuclide.clear()
        self.dlg.comboBox_nuclide.addItems(self.nuclide_list)

        # get the selected shp basin layer, get its fields and add them to the desired combobox list
        selectedLayerIndex_shp_basin = self.dlg.comboBox_shp_basin.currentIndex()
        print(f"L'indice sélectionné est : {selectedLayerIndex_shp_basin}")
        selectedLayer_shp_basin = self.layers[selectedLayerIndex_shp_basin]
        self.dlg.comboBox_ConcField.clear()
        try : # use of try method to avoid error when a DEM is wrongly selected as a the basin shp
            fields = selectedLayer_shp_basin.pendingFields()
            field_names = [field.name() for field in fields]
            self.dlg.comboBox_ConcField.addItems(field_names)
        except:
            pass


        # For denudation tab,

        self.dlg.comboBox_shp_erosion.clear()
        self.dlg.comboBox_shp_erosion.addItems(self.layer_list)
        self.dlg.comboBox_shp_erosion.currentIndexChanged[str].connect(self.shp_erosion_changed)

        # get the selected shp basin layer, get its fields and add them to the desired combobox list
        selectedLayerIndex_shp_erosion = self.dlg.comboBox_shp_erosion.currentIndex()
        selectedLayer_shp_erosion = self.layers[selectedLayerIndex_shp_erosion]
        self.dlg.comboBox_UncConc.clear()
        self.dlg.comboBox_ConcField_erosion.clear()

        try:# use of try method to avoid error when a DEM is wrongly selected as a the basin shp
            fields = selectedLayer_shp_erosion.pendingFields()
            field_names = [field.name() for field in fields]
            self.dlg.comboBox_ConcField_erosion.addItems(field_names)
            self.dlg.comboBox_UncConc.addItems(field_names)
        except:
            pass


        # show the dialog window
        self.dlg.show()

        # Run the dialog event loop
        #result = self.dlg.exec_()
        result = 0
        # See if OK was pressed
        if result:
            # architecture changed, results doesn't do anything
            pass

    def select_output_DEM_basin_file(self):
        #open the save file dialog and save the results in the proper line box
        filename, _ = QFileDialog.getSaveFileName(self.dlg, "Select output file ", "", '*.tif')
        self.dlg.lineEdit_output_DEM_basin.setText(filename)

    def select_output_qz_polygons(self):
        # open the save file dialog and save the results in the proper line box
        filename, _ = QFileDialog.getSaveFileName(self.dlg, "Select output file ", "", '*.shp')
        self.dlg.lineEdit_output_qz_polygons.setText(filename)


    def select_output_glims_polygons(self):
        # open the save file dialog and save the results in the proper line box
        filename, _ = QFileDialog.getSaveFileName(self.dlg, "Select output file ", "", '*.shp')
        self.dlg.lineEdit_output_glims_polygons.setText(filename)

    def sql_expression(self):

        try:# use of try method to avoid error when a DEM is wrongly selected as a the basin shp
            # define the expression dialog window

            # get the selected shp basin layer, get its fields and add them to the desired combobox list
            selectedLayerIndex_shp_geol = self.dlg.comboBox_shp_geol.currentIndex()
            selectedLayer_shp_geol = self.layers[selectedLayerIndex_shp_geol]
            self.dlgsql.comboBox_field.clear()
            self.dlgsql.comboBox_equation.clear()
            fields = selectedLayer_shp_geol.pendingFields()
            field_names = [field.name() for field in fields]

            # list of operator used
            list_symbols = ['=','!=','>','>=','<','<=','start with','end with']
            # add the list of operator
            self.dlgsql.comboBox_field.addItems(field_names)
            self.dlgsql.comboBox_equation.addItems(list_symbols)
            self.dlgsql.show() # show the dialog
        except:
            pass

    def check_geol(self):
        # perform the change of the desired boolean value for the corrections
        if self.dlg.checkBox_use_geol_mask.isChecked():
            self.use_geol_mask = True
        else:
            self.use_geol_mask = False

    def check_ts(self):
        # perform the change of the desired boolean value for the corrections
        if self.dlg.checkBox_use_ts_mask.isChecked():
            self.use_ts_mask = True
        else:
            self.use_ts_mask = False

    def check_VDM(self):
        # perform the change of the desired boolean value for the corrections
        if self.dlg.checkBox_VDM.isChecked():
            self.use_VDM = True
        else:
            self.use_VDM = False

    def check_glims(self):
        # perform the change of the desired boolean value for the corrections
        if self.dlg.checkBox_use_glims_mask.isChecked():
            self.use_glims_mask = True
        else:
            self.use_glims_mask = False

    def check_ifVDM_used(self):
        # perform the change of the desired boolean value for the corrections
        if self.dlg.checkBox_ifVDM.isChecked():
            self.ifVDM = True
        else:
            self.ifVDM = False

    def shp_changed(self,index):
        # if the selected basin shp changes, the dialog window is updated (same code as in the run method)
        selectedLayerIndex_shp_basin = self.dlg.comboBox_shp_basin.currentIndex()
        selectedLayer_shp_basin = self.layers[selectedLayerIndex_shp_basin]

        self.dlg.comboBox_ConcField.clear()

        try :
            fields = selectedLayer_shp_basin.pendingFields()
            field_names = [field.name() for field in fields]
            self.dlg.comboBox_ConcField.addItems(field_names)
        except :
            pass

    def shp_erosion_changed(self,index):
        # if the selected basin shp changes, the dialog window is updated (same code as in the run method)
        selectedLayerIndex_shp_erosion = self.dlg.comboBox_shp_erosion.currentIndex()
        selectedLayer_shp_erosion = self.layers[selectedLayerIndex_shp_erosion]

        self.dlg.comboBox_UncConc.clear()
        self.dlg.comboBox_ConcField_erosion.clear()

        try:
            fields = selectedLayer_shp_erosion.pendingFields()
            field_names = [field.name() for field in fields]
            self.dlg.comboBox_ConcField_erosion.addItems(field_names)
            self.dlg.comboBox_UncConc.addItems(field_names)
        except:
            pass

    def pickAlgoLS(self):
        # perform the change of the desired string value for the algorithm used in the calculation
        self.dlg.radioButton_ls.setChecked(True)
        self.dlg.radioButton_lsd.setChecked(False)
        self.dlg.radioButton_lsdp.setChecked(False)
        self.lsORlsd = 'ls'

    def pickAlgoLSD(self):
        # perform the change of the desired string value for the algorithm used in the calculation
        self.dlg.radioButton_lsd.setChecked(True)
        self.dlg.radioButton_ls.setChecked(False)
        self.dlg.radioButton_lsdp.setChecked(False)
        self.lsORlsd = 'lsd'
        
    def pickAlgoLSDp(self):
        # perform the change of the desired string value for the algorithm used in the calculation
        self.dlg.radioButton_lsdp.setChecked(True)
        self.dlg.radioButton_ls.setChecked(False)
        self.dlg.radioButton_lsd.setChecked(False)
        self.lsORlsd = 'lsdp'

    def cancel(self):
        self.dlg.close()# close the window

    def help(self):
        # open the help window
        self.dlghelp.show()

    def ok_P(self):
        # when OK on the denudation tab is pressed
        # get the content of the comboboxes and the other boxes
        output_DEM_basin = self.dlg.lineEdit_output_DEM_basin.text()

        output_qz_polygons = self.dlg.lineEdit_output_qz_polygons.text()

        express_qz_extract = self.dlg.lineEdit_express_qz_extract.text()
        # process the expression to extract the value, the operator and the field
        express_qz_extract = express_qz_extract.replace(' ','')
        if express_qz_extract == '' :
            express_qz_extract = '\"quartz\"=\'yes\''

        # get the name for the output dem, it will be used for the temp directory
        output_glims_polygons = self.dlg.lineEdit_output_glims_polygons.text()
        name = output_DEM_basin[:-4].encode('ascii','ignore')
        name = name.decode('utf-8')  # MZ Decode bytes to string
        name = name.split('/')

        # get the content of the comboboxes
        selectedLayerIndex_shp_basin = self.dlg.comboBox_shp_basin.currentIndex()
        selectedLayer_shp_basin = self.layers[selectedLayerIndex_shp_basin]

        selectedLayerIndex_DEM_raster = self.dlg.comboBox_DEM_raster.currentIndex()
        selectedLayer_DEM_raster = self.layers[selectedLayerIndex_DEM_raster]

        selectedLayerIndex_shp_geol = self.dlg.comboBox_shp_geol.currentIndex()
        selectedLayer_shp_geol = self.layers[selectedLayerIndex_shp_geol]

        selectedLayerIndex_ts_raster = self.dlg.comboBox_ts_raster.currentIndex()
        selectedLayer_ts_raster = self.layers[selectedLayerIndex_ts_raster]

        #selectedVDMIndex = self.dlg.comboBox_VDM.currentIndex()
        #selectedVDM = self.VDM_list[selectedVDMIndex]

        selectednuclideIndex = self.dlg.comboBox_nuclide.currentIndex()
        selectednuclide = self.nuclide_list[selectednuclideIndex]

        selectedLayerIndex_shp_glims = self.dlg.comboBox_shp_glims.currentIndex()
        selectedLayer_shp_glims = self.layers[selectedLayerIndex_shp_glims]

        # fields = selectedLayer_shp_basin.pendingFields()
        fields = selectedLayer_shp_basin.fields() ##MZ
        selectedConcFieldIndex = self.dlg.comboBox_ConcField.currentIndex()
        selectedConcField = fields[selectedConcFieldIndex].name()

        # get the time, used to make the temp directory
        time = datetime.datetime.now()
        time = time.isoformat()
        time = time.encode('ascii','ignore')
        if isinstance(time, bytes):  # Check if `time` is a bytes object
            time = time.decode('utf-8')  # Decode bytes to string
        time = time.replace(':','')
        time = time.replace('-', '')
        time = time.replace('T', '_')
        time = time[:15]
        temp_dir = os.path.dirname(os.path.realpath(output_DEM_basin)) + "/__temp__" + name[-1] + time

        # create a directory for the process, uses the name of the output DEM and the time
        os.mkdir(temp_dir,)

        # run pr_process , the main calculation process with the desired parameters
        pr_processes(selectedLayer_shp_basin,selectedLayer_DEM_raster,output_DEM_basin,temp_dir,selectednuclide,self.use_geol_mask,selectedLayer_shp_geol,express_qz_extract,output_qz_polygons,self.use_ts_mask,selectedLayer_ts_raster,self.use_VDM,selectedConcField,self.use_glims_mask,selectedLayer_shp_glims,output_glims_polygons,self.lsORlsd)

        # Process successful
        # self.iface.messageBar().pushMessage("Process successful", "You will find the results in the specified files; a working directory was created, you will be able to use it or delete when closing QGIS", level=QgsMessageBar.INFO, duration=20)
        self.iface.messageBar().pushMessage(
            "Process successful",
            "You will find the results in the specified files; a working directory was created, you will be able to use it or delete when closing QGIS",
            level=Qgis.Info,
            duration=20
        )
        self.dlg.close() # close the dialog window

    def ok_D(self):
        # when OK on the denudation tab is pressed
        # get the content of the combobox and the other boxes
        density = self.dlg.lineEdit_density.text()
        density = float(density)

        selectedLayerIndex_shp_erosion = self.dlg.comboBox_shp_erosion.currentIndex()
        selectedLayer_shp_erosion = self.layers[selectedLayerIndex_shp_erosion]

        fields = selectedLayer_shp_erosion.pendingFields()
        selectedConcFieldErosionIndex = self.dlg.comboBox_ConcField_erosion.currentIndex()
        selectedConcFieldErosion = fields[selectedConcFieldErosionIndex].name()

        selectedUncConcFieldIndex = self.dlg.comboBox_UncConc.currentIndex()
        selectedUncConcField = fields[selectedUncConcFieldIndex].name()


        # run denudation , the main calculation process with the desired parameters
        denudation(selectedLayer_shp_erosion,density,selectedConcFieldErosion,selectedUncConcField,self.ifVDM)


        # Process successful
        self.iface.messageBar().pushMessage("Process successful", "You will find the results in the specified file", level=QgsMessageBar.INFO, duration=20)

        self.dlg.close() # close the dialog

    def ok_sql(self):
        # when OK on the expression window is pressed
        # save the content of the comboBox and the line box to a string expression
        # The expression is then saved in the line box on the main window
        value_extract = self.dlgsql.lineEdit_value.text()

        layers = QgsProject.instance().mapLayers().values()
        selectedLayerIndex_shp_geol = self.dlg.comboBox_shp_geol.currentIndex()
        selectedLayer_shp_geol = layers[selectedLayerIndex_shp_geol]
        fields = selectedLayer_shp_geol.pendingFields()
        field_names = [field.name() for field in fields]

        list_symbols = ['=', '!=', '>', '>=', '<', '<=', 'START_With', 'END_With']

        selected_field_idx = self.dlgsql.comboBox_field.currentIndex()
        selected_eq_idx = self.dlgsql.comboBox_equation.currentIndex()
        selected_field = field_names[selected_field_idx]
        selected_eq = list_symbols[selected_eq_idx]

        sql_express = '\"' + selected_field + '\" ' + selected_eq + ' \'' + value_extract + '\''


        self.dlg.lineEdit_express_qz_extract.setText(sql_express)
        # close the dialog
        self.dlgsql.close()
        # Process successful
        self.iface.messageBar().pushMessage("Process successful", "You will find the results in the specified file", level=QgsMessageBar.INFO, duration=20)

    def cancel_sql(self):
        self.dlgsql.close() # close the window


