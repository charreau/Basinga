# -*- coding: utf-8 -*-
"""
/***************************************************************************
 BasingaDialog
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
# This creates 3 GUI dialog objects for Basinga, with the loaded .ui files

import os
from PyQt5 import QtGui,QtWidgets, uic
from qgis.utils import iface
from qgis.core import Qgis

if not os.path.exists(os.path.join( os.path.dirname(__file__), "P_E_process_dialog_base.ui")):
    test = os.path.join(
    os.path.dirname(__file__), 'P_E_process_dialog_base.ui')
    iface.messageBar().pushMessage("Problem", "No ui {}".format(test), level=Qgis.Warning, duration=10)

    
FORM_CLASS, _ = uic.loadUiType(os.path.join( os.path.dirname(__file__), "P_E_process_dialog_base.ui"))
class BasingaDialog(QtWidgets.QDialog, FORM_CLASS):
    def __init__(self, parent=None):
        """Constructor."""
        super(BasingaDialog, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)


FORM_CLASS_sql, _sql = uic.loadUiType(os.path.join( os.path.dirname(__file__), "P_E_process_dialog_sql.ui"))
class SQLDialog(QtWidgets.QDialog, FORM_CLASS_sql):
    def __init__(self, parent=None):
        """Constructor."""
        super(SQLDialog, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)


FORM_CLASS_help, _help = uic.loadUiType(os.path.join( os.path.dirname(__file__), "P_E_process_dialog_help.ui"))
class HelpDialog(QtWidgets.QDialog, FORM_CLASS_help):
    def __init__(self, parent=None):
        """Constructor."""
        super(HelpDialog, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)