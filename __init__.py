version="1.2.1"
import sys
path=__path__[0]+'/'
sys.path.append(path)

import simpleLoader
simpleLoader.path=path
from simpleLoader import startMC
