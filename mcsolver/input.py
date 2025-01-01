import numpy as np
import Lattice as lat
import mcMain as mc
import matplotlib.pyplot as plt
import time

## magnetic crystal part
LMatrix=[[1,0,0],
         [0,1,0],
         [0,0,1]]
# magnetic orbitals in fractional coordinates
pos=[[0,0,0]] 
# spin number
Spin=[1]
# single ion anisotropy
D=np.array([[0.0,0.0,-0.0]])
# couplings   #source #target #edge  #J(meV) negative for FM coupling
bond1=lat.Bond(0,0,np.array([1,0,0]),-1,-1,-1, 0, 0, 0, 0, 0, 0, True)
bond2=lat.Bond(0,0,np.array([0,1,0]),-1,-1,-1, 0, 0, 0, 0, 0, 0, True)

bondList=[bond1,bond2]

time0=time.time()
mcslave=mc.MC(0,LMatrix,pos=pos,S=Spin,D=D,bondList=bondList,
              Lx=16,Ly=16,Lz=1,
              T=0.2,ki_s=0,ki_t=0,ki_overLat=[2,0,0],
              orbGroupList=[[0]],groupInSC=True,h=0.0,dipoleAlpha=0.0,On=1,spinFrame=0)
# Check MC part
data=mcslave.mainLoopViaCLib_On(nsweep=80000,nthermal=40000,ninterval=1,algo='Wolff',On=3,flunc=0.0)


# Spin-wave part

print('time elapsed: %.3f s'%(time.time()-time0))

