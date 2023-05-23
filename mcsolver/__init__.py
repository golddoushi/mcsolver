__version__="3.0.0"
from . import win
import os

libpath=os.path.join(os.path.dirname(__file__),'lib')
for lib in os.listdir(libpath):
    if 'ising' in lib:
        isingLibPath=os.path.join(libpath,lib)
    if 'xy' in lib:
        xyLibPath=os.path.join(libpath,lib)
    if 'heisenberg' in lib:
        heisenbergLibPath=os.path.join(libpath,lib)

win.libPool=[isingLibPath,xyLibPath,heisenbergLibPath]

def loadMC(rpath): # interface to core codes avoiding GUI
    win.startSimulation(updateGUI=False,rpath=rpath)

