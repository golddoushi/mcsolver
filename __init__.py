from tkinter import Tk
from multiprocessing import Pool, freeze_support
import numpy as np
import time
import guiMain as gui
import Lattice as lat
import mcMain as mc

global bondList,LMatrix,pos,nsweep,nthermal,Lx,Ly,Lz,algorithm

def startMC(param): # start MC for Ising model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,nsweep,nthermal,Lx,Ly,Lz,algorithm=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz)
    mData, eData=np.array(mcslave.mainLoopViaCLib(nsweep=nsweep,nthermal=nthermal,algo=algorithm))
    mData=abs(mData)/Lx/Ly/Lz
    eData/=(Lx*Ly*Lz)
    return ID, T, mData, eData

def startMCForOn(param): # start MC for O(n) model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,nsweep,nthermal,Lx,Ly,Lz,algorithm,On=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz)
    mData, eData=mcslave.mainLoopViaCLib_On(nsweep=nsweep,nthermal=nthermal,algo=algorithm,On=On)
    mData=abs(mData)/Lx/Ly/Lz
    eData/=(Lx*Ly*Lz)
    return ID, T, mData, eData

def startSimulaton():
    time0=time.time()
    global bondList,LMatrix,pos,nsweep,nthermal,Lx,Ly,Lz,algorithm
    gui.submitBtn.config(state='disabled')
    # get lattice
    a1=gui.latticeGui[0].report()
    a2=gui.latticeGui[1].report()
    a3=gui.latticeGui[2].report()
    LMatrix=np.array([a1,a2,a3])
    print('Lattice matrix:')
    print(LMatrix)

    # get supercell size
    Lx, Ly, Lz=[int(x) for x in gui.supercellGui.report()]
    print('supercell:')
    print(Lx,Ly,Lz)

    # get oribtal position and spin state and onsite-anisotropy
    pos=np.array([ele[3] for ele in gui.OrbListBox.infoData])
    S=[ele[2] for ele in gui.OrbListBox.infoData]
    DList=[ele[4] for ele in gui.OrbListBox.infoData]
    for ipos, iS, iD in zip(pos,S,DList):
        print('positions:',ipos,'Spin:',iS,'onsite-Anisotropy:',iD)

    # get bonds
    bondList=[lat.Bond(bond_data[2][0],bond_data[2][1],\
                       np.array([int(x) for x in bond_data[2][2]]),\
                       bond_data[1][0],bond_data[1][1],bond_data[1][2]) \
                        for bond_data in gui.BondBox.infoData]
        
    print('bonds:')
    print(bondList)

    # get TList
    T0, T1, nT=gui.TListGui.report()
    TList=np.linspace(T0,T1,int(nT))
    print('Temperature:')
    print(TList)

    # get thermalizations and sweeps
    nthermal, nsweep = [int(x) for x in gui.MCparamGui.report()]
    print('thermalizations and sweeps:')
    print(nthermal, nsweep)

    # get model and algorithm
    modelType = gui.modelGui.get()
    print('Model:',modelType)
    algorithm = gui.algorithmGui.get()
    print('Algorithm:',algorithm)

    # get ncores
    ncores= int(gui.coreGui.report()[0])
    print('using %d cores'%ncores)

    # model and algorithm branches
    # Ising model settings
    if(modelType=='Ising'):
        if algorithm!='Metroplis' and algorithm!='Wolff':
            print('For now, only Metroplis and Wolff algorithm is supported for Ising model')
            gui.submitBtn.config(state='normal')
            return
        
        paramPack=[]
        for iT, T in enumerate(TList):
            paramPack.append([iT,T,bondList,LMatrix,pos,S,DList,nsweep,nthermal,Lx,Ly,Lz,algorithm])
        
        TResult=[];magResult=[];susResult=[];energyResult=[];capaResult=[]
        pool=Pool(processes=ncores)
        for result in pool.imap_unordered(startMC,paramPack):
            ID, T, mData, eData =result
            TResult.append(T)
            magResult.append(np.mean(mData))
            susResult.append(np.std(mData))
            energyResult.append(np.mean(eData))
            capaResult.append(np.std(eData))
        pool.close()
        gui.updateResultViewer(TList=TResult, magList=magResult, susList=susResult)
    # continuous model settings
    elif(modelType=='XY' or modelType=='Heisenberg'):
        for bond in bondList:
            bond.On=True # switch on the vector type bonding
        if algorithm!='Metroplis' and algorithm!='Wolff':
            print('For now, only Metroplis and Wolff algorithm is supported for O(n) model')
            gui.submitBtn.config(state='normal')
            return

        On=2 if modelType=='XY' else 3
        paramPack=[]
        for iT, T in enumerate(TList):
            paramPack.append([iT,T,bondList,LMatrix,pos,S,DList,nsweep,nthermal,Lx,Ly,Lz,algorithm,On])

        TResult=[];magResult=[];susResult=[];energyResult=[];capaResult=[]
        pool=Pool(processes=ncores)
        for result in pool.imap_unordered(startMCForOn,paramPack):
            ID, T, mData, eData =result
            TResult.append(T)
            magResult.append(np.mean(mData))
            susResult.append(np.std(mData))
            energyResult.append(np.mean(eData))
            capaResult.append(np.std(eData))
        pool.close()
        gui.updateResultViewer(TList=TResult, magList=magResult, susList=susResult)
    else:
        print("for now only Ising, XY and Heisenberg model is supported")
        gui.submitBtn.config(state='normal')
        return

    # writting result file
    f=open('./result.txt','w')
    f.write('#Temp #Spin    #Susc      #energy  #capacity\n')
    for T, mag, sus, energy, capa in zip(TResult, magResult, susResult, energyResult, capaResult):
        f.write('%.3f %.6f %.6f %.6f %.6f\n'%(T, mag, sus, energy, capa))
    f.close()
    gui.submitBtn.config(state='normal')
    print("time elapsed %.3f s"%(time.time()-time0))
    return
            
if __name__ == '__main__': # crucial for multiprocessing in Windows
    freeze_support()
    app=Tk(className='mc solver v1.0')
    gui.loadEverything(app,startSimulaton)
    app.mainloop()