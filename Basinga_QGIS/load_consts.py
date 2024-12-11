import numpy as np
import scipy.io
import scipy.interpolate
import sys, os
import json

path = os.path.dirname(sys.argv[0])+ '/consts_py.mat'
print path
mat = scipy.io.loadmat(path)

# Open, load, save and return consts_py.mat in a dictionary
path = 'C:\Users\CHARREAU\.qgis2\python\plugins\Basinga\consts_py.mat'
mat = scipy.io.loadmat(path)
struc = mat['consts']

listVar = struc.dtype.names
valGlob = struc[0, 0]
constsDict = {}
for name in listVar:
    constsDict[name] = valGlob[name]

#np.save('C:\Users\CHARREAU\.qgis2\python\plugins\Basinga\consts', constsDict) 
file=open('C:\Users\CHARREAU\.qgis2\python\plugins\Basinga\consts.py','w')

file.write(str(constsDict))

file.flush()

