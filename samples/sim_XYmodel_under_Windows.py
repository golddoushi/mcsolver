from multiprocessing import freeze_support
import mcsolver

if __name__ == '__main__': # crucial for multiprocessing in Windows
    freeze_support()
    mcsolver.loadMC("./Square_XY_isotropic")