from tkinter import Tk
from multiprocessing import Pool, freeze_support
import numpy as np
import time
import guiMain as gui
import Lattice as lat
import mcMain as mc
import fileio as io

global path
path='./'

def startMC(param): # start MC for Ising model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,h,nsweep,nthermal,ninterval,Lx,Ly,Lz,algorithm,GcOrb,dipoleAlpha=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz,ki_s=GcOrb[0][0],ki_t=GcOrb[0][1],ki_overLat=GcOrb[1],h=h,dipoleAlpha=dipoleAlpha)
    spin_i, spin_j, spin_ij, autoCorr, E, E2, U4=mcslave.mainLoopViaCLib(nsweep=nsweep,nthermal=nthermal,ninterval=ninterval,algo=algorithm)
    #mData=abs(mData)/Lx/Ly/Lz
    #eData/=(Lx*Ly*Lz)
    #print("<ij>=",np.mean(corr))
    return ID, T, h, spin_i, spin_j, spin_ij, autoCorr, E, E2, U4, mcslave.totOrbs

def startMCForOn(param): # start MC for O(n) model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,h,nsweep,nthermal,ninterval,Lx,Ly,Lz,algorithm,On,GcOrb,dipoleAlpha=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz,ki_s=GcOrb[0][0],ki_t=GcOrb[0][1],ki_overLat=GcOrb[1],h=h,dipoleAlpha=dipoleAlpha,On=On)
    spin_i, spin_j, spin_ij, autoCorr, E, E2, U4=mcslave.mainLoopViaCLib_On(nsweep=nsweep,nthermal=nthermal,ninterval=ninterval,algo=algorithm,On=On)
    #mData=abs(mData)/Lx/Ly/Lz
    #eData/=(Lx*Ly*Lz)
    return ID, T, h, spin_i, spin_j, spin_ij, autoCorr, E, E2, U4, mcslave.totOrbs

def startSimulation(updateGUI=True, rpath=''):
    time0=time.time()
    with open('./out','w') as fout:
        fout.write('#T #H\n')

    if updateGUI:
        gui.submitBtn.config(state='disabled')
        io.collectParam()
    else:
        io.loadParam(updateGUI=False, rpath=rpath)

    TList=np.linspace(io.T0,io.T1,io.nT)
    HList=np.linspace(io.H0,io.H1,io.nH)
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
            if updateGUI: gui.submitBtn.config(state='normal')
            return
        
        paramPack=[]
        for iH, H in enumerate(HList):
            for iT, T in enumerate(TList):
                paramPack.append([iH*len(TList)+iT,T,bondList,LMatrix,pos,io.S,io.DList,H,io.nsweep,io.nthermal,io.ninterval,io.LPack[0],io.LPack[1],io.LPack[2],io.algorithm,
                                 io.GcOrb,io.dipoleAlpha])
        
        TResult=[];HResult=[];SpinIResult=[];SpinJResult=[];susResult=[];energyResult=[];capaResult=[];u4Result=[];autoCorrResult=[]
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
                ID, T, h, spin_i, spin_j, spin_ij, autoCorr, E, E2, U4, N=result
                TResult.append(T)
                HResult.append(h)
                SpinIResult.append(spin_i)
                SpinJResult.append(spin_j)
                susResult.append((spin_ij-spin_i*spin_j)/T)
                autoCorrResult.append(autoCorr)
                energyResult.append(E)
                capaResult.append((E2-E*E)/T**2)
                u4Result.append(U4)
            pool.close()
        if updateGUI: gui.updateResultViewer(TList=HResult, magList=[(si+sj)/2 for si,sj in zip(SpinIResult,SpinJResult)], susList=capaResult)
    # continuous model settings
    elif(io.modelType=='XY' or io.modelType=='Heisenberg'):
        for bond in bondList:
            bond.On=True # switch on the vector type bonding
        if io.algorithm!='Metroplis' and io.algorithm!='Wolff':
            print('For now, only Metroplis and Wolff algorithm is supported for O(n) model')
            if updateGUI: gui.submitBtn.config(state='normal')
            return

        On=2 if io.modelType=='XY' else 3
        paramPack=[]
        for iH, H in enumerate(HList):
            for iT, T in enumerate(TList):
                paramPack.append([iH*len(TList)+iT,T,bondList,LMatrix,pos,io.S,io.DList,H,io.nsweep,io.nthermal,io.ninterval,io.LPack[0],io.LPack[1],io.LPack[2],io.algorithm,On,
                                  io.GcOrb,io.dipoleAlpha])

        TResult=[];HResult=[];SpinIResult=[];SpinJResult=[];susResult=[];energyResult=[];capaResult=[];u4Result=[];autoCorrResult=[]
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
                ID, T, h, spin_i, spin_j, spin_ij, autoCorr, E, E2, U4, N =result
                TResult.append(T)
                HResult.append(h)
                SpinIResult.append(np.sqrt(sum(spin_i*spin_i)))
                SpinJResult.append(np.sqrt(sum(spin_j*spin_j)))
                susResult.append((spin_ij-np.dot(spin_i,spin_j))/T)
                autoCorrResult.append(autoCorr)
                energyResult.append(E)
                capaResult.append((E2-E*E)/T**2)
                u4Result.append(U4)
            pool.close()
        if updateGUI: gui.updateResultViewer(TList=HResult, magList=SpinIResult, susList=capaResult)
    else:
        print("for now only Ising, XY and Heisenberg model is supported")
        if updateGUI: gui.submitBtn.config(state='normal')
        return

    # writting result file
    f=open('./result.txt','w')
    f.write('#Temp #<Si>    #<Sj>    #Susc    #energy   #capacity #Binder cumulante #auto-corr.\n')
    for T, si, sj, sus, energy, capa, u4, autoCorr in zip(TResult, SpinIResult, SpinJResult, susResult, energyResult, capaResult, u4Result, autoCorrResult):
        f.write('%.3f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n'%(T, si, sj, sus, energy, capa, u4, autoCorr))
    f.close()
    if updateGUI: gui.submitBtn.config(state='normal')
    print("time elapsed %.3f s"%(time.time()-time0))
    return
            
if __name__ == '__main__': # crucial for multiprocessing in Windows
    freeze_support()
    app=Tk(className='mc solver v1.0')
    gui.loadEverything(app,startSimulation)
    app.mainloop()