# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Basinga
                                 A QGIS plugin
 Caculate 10Be production rates and erosion
                             -------------------
        begin                : 2016-06-21
        copyright            : (C) 2016 by CRPG - Nancy France
        email                : charreau@crpg.cnrs-nancy.fr
        git sha              :
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load Basinga class from file Basinga.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """

    from .P_E_process import Basinga
    return Basinga(iface)
