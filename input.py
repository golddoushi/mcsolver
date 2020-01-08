import numpy as np
import Lattice as lat
import mcMain as mc
import matplotlib.pyplot as plt
import time
## magnetic crystal part
LMatrix=[[1,0,0],
         [-0.5,0.8660254,0],
         [0,0,1]]
# magnetic orbitals in fractional coordinates
pos=[[2./3,1./3,0],
     [1./3,2./3,0]] 
# spin number
Spin=[1.5,1.5]
# single ion anisotropy
D=[[-2.4144,0.0,0.0],
   [-2.4144,0.0,0.0]]
# couplings   #source #target #edge  #J(meV) negative for FM coupling
bond1=lat.Bond(0,1,np.array([0,0,0]),-31.3542,-29.857,-29.857, False)
bond2=lat.Bond(0,1,np.array([1,0,0]),-31.3542,-29.857,-29.857, False)
bond3=lat.Bond(0,1,np.array([0,-1,0]),-31.3542,-29.857,-29.857, False)

bondList=[bond1,bond2,bond3
         ]


time0=time.time()
mcslave=mc.MC(0,LMatrix,pos=pos,S=Spin,D=D,bondList=bondList,T=50,Lx=16,Ly=16,Lz=1)
#mcslave.mainLoopViaCLib_On(nsweep=1,nthermal=1,algo='Wolff',On=3,flunc=0.0)
mcslave.mainLoopViaCLib(nsweep=1,nthermal=1,algo='Wolff')
print('time elapsed: %.3f s'%(time.time()-time0))

'''
TList=np.linspace(45,55,10)
U4_list=[]
for L in [16,32]:
    U4_L=[]
    for T in TList:
        mcslave=mc.MC(0,LMatrix,pos=pos,S=Spin,D=D,bondList=bondList,T=T,Lx=L,Ly=L,Lz=1)
        mData, eData=np.array(mcslave.mainLoopViaCLib_On(nsweep=320000,nthermal=80000,algo='Wolff',On=3,flunc=0.0))
        U4_L.append(np.mean(mData*mData)**2/np.mean(mData**4))
        #U4_L.append(np.mean(mData))
    U4_list.append(U4_L)
    
plt.scatter(TList,U4_list[0])
plt.scatter(TList,U4_list[1])
plt.show()
'''