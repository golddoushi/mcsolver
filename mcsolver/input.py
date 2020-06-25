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
D=np.array([[0.0,0.0,0.0]])
# couplings   #source #target #edge  #J(meV) negative for FM coupling
bond1=lat.Bond(0,0,np.array([1,0,0]),-1,-1,-1, False)
bond2=lat.Bond(0,0,np.array([0,1,0]),-1,-1,-1, False)

bondList=[bond1,bond2]

time0=time.time()
mcslave=mc.MC(0,LMatrix,pos=pos,S=Spin,D=D,bondList=bondList,T=2,Lx=32,Ly=32,Lz=1,ki_s=0,ki_t=0,ki_overLat=[0,0,0],h=0.01,dipoleAlpha=0.0,On=1,spinFrame=1)
#data=mcslave.mainLoopViaCLib_On(nsweep=360000,nthermal=80000,ninterval=1,algo='Wolff',On=3,flunc=0.0,binGraph=False)
data=mcslave.mainLoopViaCLib(nsweep=20000,nthermal=10000,ninterval=0,algo='Metroplis')
print(data)
#mean=np.mean(abs(np.array(totSpin)))/mcslave.totOrbs
#print(np.mean(abs(np.array(corr)))-mean**2)
print('time elapsed: %.3f s'%(time.time()-time0))

'''
TList=np.linspace(2.1,2.4,11)
U4_list=[]
for L in [16,32]:
    U4_L=[]
    for T in TList:
        mcslave=mcslave=mc.MC(0,LMatrix,pos=pos,S=Spin,D=D,bondList=bondList,T=T,Lx=L,Ly=L,Lz=1,ki_s=0,ki_t=0,ki_overLat=[0,0,0],h=0.0)
        totSpin,energy,corr=mcslave.mainLoopViaCLib(nsweep=1000,nthermal=500,ninterval=100,algo='Wolff')
        mean=np.mean(abs(np.array(totSpin)))/mcslave.totOrbs
        U4_L.append((np.mean(abs(np.array(corr)))-mean**2)/T)
        #U4_L.append(np.std(abs(np.array(totSpin))/mcslave.totOrbs)**2/T)
        #U4_L.append(np.mean(mData))
    U4_list.append(U4_L)
    
plt.scatter(TList,U4_list[0])
plt.scatter(TList,U4_list[1])
plt.show()
'''
