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
Spin=[1.]
D=[[0.0,0.0,0.0]]  # single ion anisotropy
# couplings   #source #target #edge  #J(meV) negative for FM coupling
bond1=lat.Bond(0,0,np.array([1,0,0]),-13.,-10.0,-10.0, True)
bond2=lat.Bond(0,0,np.array([0,1,0]),-13.,-10.0,-10.0, True)
bondList=[bond1,bond2]


#time0=time.time()
#mcslave=mc.MC(0,LMatrix,pos=pos,S=Spin,D=D,bondList=bondList,T=9.0,Lx=32,Ly=32,Lz=1)
#mcslave.mainLoopViaCLib_On(nsweep=40000,nthermal=20000,algo='Wolff',On=3,flunc=0.0)
#print('time elapsed: %.3f s'%(time.time()-time0))

'''
TList=np.linspace(9,11,10)
U4_list=[]
for L in [32,64]:
    U4_L=[]
    for T in TList:
        mcslave=mc.MC(0,LMatrix,pos=pos,S=Spin,D=D,bondList=bondList,T=T,Lx=L,Ly=L,Lz=1)
        mData, eData=np.array(mcslave.mainLoopViaCLib_On(nsweep=40000,nthermal=20000,algo='Wolff',On=3,flunc=0.0))
        U4_L.append(np.mean(mData*mData)**2/np.mean(mData**4))
    U4_list.append(U4_L)
    
plt.plot(TList,U4_list[0])
plt.plot(TList,U4_list[1])
plt.show()
'''