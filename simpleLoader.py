import fileio as io
from multiprocessing import Pool, freeze_support
import numpy as np
import time
import Lattice as lat
import mcMain as mc
import fileio as io

global path
path='./'

def startMCForIsing(param): # start MC for Ising model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,nsweep,nthermal,ninterval,Lx,Ly,Lz,algorithm,GcOrb=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz,ki_s=GcOrb[0][0],ki_t=GcOrb[0][1],ki_overLat=GcOrb[1])
    spin_i, spin_j, spin_ij, E, E2, U4=mcslave.mainLoopViaCLib(nsweep=nsweep,nthermal=nthermal,ninterval=ninterval,algo=algorithm)
    #mData=abs(mData)/Lx/Ly/Lz
    #eData/=(Lx*Ly*Lz)
    #print("<ij>=",np.mean(corr))
    return ID, T, spin_i, spin_j, spin_ij, E, E2, U4, mcslave.totOrbs

def startMCForOn(param): # start MC for O(n) model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,nsweep,nthermal,ninterval,Lx,Ly,Lz,algorithm,On,GcOrb=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz,ki_s=GcOrb[0][0],ki_t=GcOrb[0][1],ki_overLat=GcOrb[1],h=0.1)
    spin_i, spin_j, spin_ij, E, E2, U4=mcslave.mainLoopViaCLib_On(nsweep=nsweep,nthermal=nthermal,ninterval=ninterval,algo=algorithm,On=On)
    #mData=abs(mData)/Lx/Ly/Lz
    #eData/=(Lx*Ly*Lz)
    return ID, T, spin_i, spin_j, spin_ij, E, E2, U4, mcslave.totOrbs

def startMC(paramPath):
    io.loadParam(updateGUI=False, rpath=paramPath)

    time0=time.time()
    
    TList=np.linspace(io.T0,io.T1,io.nT)
    bondList=[lat.Bond(bond_data[0],bond_data[1],                 # source and target  
                       np.array([int(x) for x in bond_data[2]]),  # over lat.
                       bond_data[3],bond_data[4],bond_data[5],    # strength
                       True if io.modelType!='Ising' else False)  # On
                        for bond_data in io.bondList]             # ergodic
    LMatrix=np.array(io.LMatrix)
    pos=np.array(io.pos)
    
    if(io.modelType=='Ising'):
        if io.algorithm!='Metroplis' and io.algorithm!='Wolff':
            print('For now, only Metroplis and Wolff algorithm is supported for Ising model')
            return
        
        paramPack=[]
        for iT, T in enumerate(TList):
            paramPack.append([iT,T,bondList,LMatrix,pos,io.S,io.DList,io.nsweep,io.nthermal,io.ninterval,io.LPack[0],io.LPack[1],io.LPack[2],io.algorithm,
                              io.GcOrb])
        
        TResult=[];magResult=[];susResult=[];energyResult=[];capaResult=[];u4Result=[]
        while(True): # using pump strategy to reduce the costs of RAM
            if len(paramPack)==0:
                break
            # pump tasks
            paramPack_tmp=[]
            for index in range(io.ncores):
                paramPack_tmp.append(paramPack.pop(0))
                if len(paramPack)==0:
                    break
            pool=Pool(processes=io.ncores)
            for result in pool.imap_unordered(startMCForIsing,paramPack_tmp):
                ID, T, spin_i, spin_j, spin_ij, E, E2, U4, N=result
                TResult.append(T)
                magResult.append(spin_i)
                susResult.append((spin_ij-spin_i*spin_j)/T)
                energyResult.append(E)
                capaResult.append((E2-E*E)/T**2)
                u4Result.append(U4)
            pool.close()
    # continuous model settings
    elif(io.modelType=='XY' or io.modelType=='Heisenberg'):
        for bond in bondList:
            bond.On=True # switch on the vector type bonding
        if io.algorithm!='Metroplis' and io.algorithm!='Wolff':
            print('For now, only Metroplis and Wolff algorithm is supported for O(n) model')
            return

        On=2 if io.modelType=='XY' else 3
        paramPack=[]
        for iT, T in enumerate(TList):
            paramPack.append([iT,T,bondList,LMatrix,pos,io.S,io.DList,io.nsweep,io.nthermal,io.ninterval,io.LPack[0],io.LPack[1],io.LPack[2],io.algorithm,On,
                              io.GcOrb])

        TResult=[];magResult=[];susResult=[];energyResult=[];capaResult=[];u4Result=[]
        while(True): # using pump strategy to reduce the costs of RAM
            if len(paramPack)==0:
                break
            # pump tasks
            paramPack_tmp=[]
            for index in range(io.ncores):
                paramPack_tmp.append(paramPack.pop(0))
                if len(paramPack)==0:
                    break
            pool=Pool(processes=io.ncores)
            for result in pool.imap_unordered(startMCForOn,paramPack_tmp):
                ID, T, spin_i, spin_j, spin_ij, E, E2, U4, N =result
                TResult.append(T)
                magResult.append(np.sqrt(sum(spin_i*spin_i)))
                susResult.append((spin_ij-np.dot(spin_i,spin_j))/T)
                energyResult.append(E)
                capaResult.append((E2-E*E)/T**2)
                u4Result.append(U4)
            pool.close()
    else:
        print("for now only Ising, XY and Heisenberg model is supported")
        return

    # writting result file
    f=open('./result.txt','w')
    f.write('#Temp #Spin    #Susc      #energy  #capacity #Binder cumulante\n')
    for T, mag, sus, energy, capa, u4 in zip(TResult, magResult, susResult, energyResult, capaResult, u4Result):
        f.write('%.3f %.6f %.6f %.6f %.6f %.6f\n'%(T, mag, sus, energy, capa, u4))
    f.close()
    print("time elapsed %.3f s"%(time.time()-time0))
    return

if __name__ == '__main__':
    freeze_support()
    startMC('./samples/1')