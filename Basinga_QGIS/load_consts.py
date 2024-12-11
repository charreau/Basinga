import numpy as np
import scipy.io
import scipy.interpolate
import sys, os
import json

path_local = os.path.dirname(sys.argv[0])
path_mat_file  = os.path.join(path_local, "/consts_py.mat")
path_dict_file = os.path.join(path_local, "/consts.py")
print(path_mat_file)
mat = scipy.io.loadmat(path_mat_file)

# Open, load, save and return consts_py.mat in a dictionary
# path = 'C:\Users\CHARREAU\.qgis2\python\plugins\Basinga\consts_py.mat'
# path='C:\Users\zakari\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins\Basinga_QGIS\consts.py.mat'
# mat = scipy.io.loadmat(path)
struc = mat['consts']

listVar = struc.dtype.names
valGlob = struc[0, 0]
constsDict = {}
for name in listVar:
    constsDict[name] = valGlob[name]

#np.save('C:\Users\CHARREAU\.qgis2\python\plugins\Basinga\consts', constsDict) 
# file=open('C:\Users\CHARREAU\.qgis2\python\plugins\Basinga\consts.py','w')
file=open(path_dict_file,'w')

file.write(str(constsDict))

file.flush()

