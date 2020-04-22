from tkinter import Tk
from multiprocessing import Pool, freeze_support
import numpy as np
import time
import guiMain as gui
import Lattice as lat
import mcMain as mc
import fileio as io

def startMC(param): # start MC for Ising model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,nsweep,nthermal,ninterval,Lx,Ly,Lz,algorithm=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz)
    spin_i, spin_j, spin_ij, E, E2=mcslave.mainLoopViaCLib(nsweep=nsweep,nthermal=nthermal,ninterval=ninterval,algo=algorithm)
    #mData=abs(mData)/Lx/Ly/Lz
    #eData/=(Lx*Ly*Lz)
    #print("<ij>=",np.mean(corr))
    return ID, T, spin_i, spin_j, spin_ij, E, E2, mcslave.totOrbs

def startMCForOn(param): # start MC for O(n) model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,nsweep,nthermal,ninterval,Lx,Ly,Lz,algorithm,On=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz,h=0.1)
    spin_i, spin_j, spin_ij, E, E2=mcslave.mainLoopViaCLib_On(nsweep=nsweep,nthermal=nthermal,ninterval=ninterval,algo=algorithm,On=On)
    #mData=abs(mData)/Lx/Ly/Lz
    #eData/=(Lx*Ly*Lz)
    return ID, T, spin_i, spin_j, spin_ij, E, E2, mcslave.totOrbs

def startSimulaton():
    time0=time.time()
    gui.submitBtn.config(state='disabled')
    io.collectParam()
    
    TList=np.linspace(io.T0,io.T1,io.nT)
    bondList=[lat.Bond(bond_data[0],bond_data[1],\
                       np.array([int(x) for x in bond_data[2]]),\
                       bond_data[3],bond_data[4],bond_data[5]) \
                        for bond_data in io.bondList]
    LMatrix=np.array(io.LMatrix)
    pos=np.array(io.pos)
    
    if(io.modelType=='Ising'):
        if io.algorithm!='Metroplis' and io.algorithm!='Wolff':
            print('For now, only Metroplis and Wolff algorithm is supported for Ising model')
            gui.submitBtn.config(state='normal')
            return
        
        paramPack=[]
        for iT, T in enumerate(TList):
            paramPack.append([iT,T,bondList,LMatrix,pos,io.S,io.DList,io.nsweep,io.nthermal,io.ninterval,io.LPack[0],io.LPack[1],io.LPack[2],io.algorithm])
        
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
            for result in pool.imap_unordered(startMC,paramPack_tmp):
                ID, T, spin_i, spin_j, spin_ij, E, E2, N=result
                TResult.append(T)
                magResult.append(spin_i)
                susResult.append((spin_ij-spin_i*spin_j)/T)
                energyResult.append(E)
                capaResult.append((E2-E*E)/T**2)
                u4Result.append(0.)
            pool.close()
        gui.updateResultViewer(TList=TResult, magList=magResult, susList=capaResult)
    # continuous model settings
    elif(io.modelType=='XY' or io.modelType=='Heisenberg'):
        for bond in bondList:
            bond.On=True # switch on the vector type bonding
        if io.algorithm!='Metroplis' and io.algorithm!='Wolff':
            print('For now, only Metroplis and Wolff algorithm is supported for O(n) model')
            gui.submitBtn.config(state='normal')
            return

        On=2 if io.modelType=='XY' else 3
        paramPack=[]
        for iT, T in enumerate(TList):
            paramPack.append([iT,T,bondList,LMatrix,pos,io.S,io.DList,io.nsweep,io.nthermal,io.ninterval,io.LPack[0],io.LPack[1],io.LPack[2],io.algorithm,On])

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
                ID, T, spin_i, spin_j, spin_ij, E, E2, N =result
                TResult.append(T)
                magResult.append(np.sqrt(sum(spin_i*spin_i)))
                susResult.append((spin_ij-np.dot(spin_i,spin_j))/T)
                energyResult.append(E)
                capaResult.append((E2-E*E)/T**2)
                u4Result.append(0.)
            pool.close()
        gui.updateResultViewer(TList=TResult, magList=magResult, susList=capaResult)
    else:
        print("for now only Ising, XY and Heisenberg model is supported")
        gui.submitBtn.config(state='normal')
        return

    # writting result file
    f=open('./result.txt','w')
    f.write('#Temp #Spin    #Susc      #energy  #capacity #Binder cumulante\n')
    for T, mag, sus, energy, capa, u4 in zip(TResult, magResult, susResult, energyResult, capaResult, u4Result):
        f.write('%.3f %.6f %.6f %.6f %.6f %.6f\n'%(T, mag, sus, energy, capa, u4))
    f.close()
    gui.submitBtn.config(state='normal')
    print("time elapsed %.3f s"%(time.time()-time0))
    return
            
if __name__ == '__main__': # crucial for multiprocessing in Windows
    freeze_support()
    app=Tk(className='mc solver v1.0')
    gui.loadEverything(app,startSimulaton)
    app.mainloop()