__version__="1.2.1"
import sys
path=__path__[0]+'/'
sys.path.append(path)
import win
win.path=path

def loadMC(rpath): # interface to core codes avoiding GUI
    win.startSimulation(updateGUI=False,rpath=rpath)

