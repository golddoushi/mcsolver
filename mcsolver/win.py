from tkinter import Tk
from multiprocessing import Pool, freeze_support
import numpy as np
import time

try:
    from . import guiMain as gui
    from . import Lattice as lat
    from . import mcMain as mc
    from . import fileio as io
except:
    import guiMain as gui
    import Lattice as lat
    import mcMain as mc
    import fileio as io

global path, settingFileVersion
libPool=[None,None,None]
settingFileVersion=2.3

def startMC(param): # start MC for Ising model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,h,nsweep,nthermal,ninterval,Lx,Ly,Lz,algorithm,GcOrb,orbGroupList,groupInSC,dipoleAlpha,spinFrame=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz,ki_s=GcOrb[0][0],ki_t=GcOrb[0][1],ki_overLat=GcOrb[1],orbGroupList=orbGroupList,groupInSC=groupInSC,h=h,dipoleAlpha=dipoleAlpha,spinFrame=spinFrame)
    spin_i, spin_j, spin_ij, autoCorr, E, E2, U4=mcslave.mainLoopViaCLib(nsweep=nsweep,nthermal=nthermal,ninterval=ninterval,algo=algorithm)
    #mData=abs(mData)/Lx/Ly/Lz
    #eData/=(Lx*Ly*Lz)
    #print("<ij>=",np.mean(corr))
    return ID, T, h, spin_i, spin_j, spin_ij, autoCorr, E, E2, U4, mcslave.totOrbs

def startMCForOn(param): # start MC for O(n) model
    # unzip all global parameters for every processing
    ID, T, bondList,LMatrix,pos,S,DList,h,nsweep,nthermal,ninterval,Lx,Ly,Lz,algorithm,On,GcOrb,orbGroupList,groupInSC,dipoleAlpha,spinFrame=param
    mcslave=mc.MC(ID,LMatrix,pos=pos,S=S,D=DList,bondList=bondList,T=T,Lx=Lx,Ly=Ly,Lz=Lz,ki_s=GcOrb[0][0],ki_t=GcOrb[0][1],ki_overLat=GcOrb[1],orbGroupList=orbGroupList,groupInSC=groupInSC,h=h,dipoleAlpha=dipoleAlpha,On=On,spinFrame=spinFrame)
    spin_i, spin_j, spin_ij, autoCorr, E, E2, U4=mcslave.mainLoopViaCLib_On(nsweep=nsweep,nthermal=nthermal,ninterval=ninterval,algo=algorithm,On=On)
    #mData=abs(mData)/Lx/Ly/Lz
    #eData/=(Lx*Ly*Lz)
    return ID, T, h, spin_i, spin_j, spin_ij, autoCorr, E, E2, U4, mcslave.totOrbs

def startSimulation(updateGUI=True, rpath=''):
    time0=time.time()
    # clean possible existed files
    with open('./out','w') as fout:
        fout.write('#T #H\n')
    with open('./spinDotSpin.txt','w') as fout:
        fout.write('#T #H\n')

    if updateGUI:
        gui.submitBtn.config(state='disabled')
        io.collectParam()
    else:
        if not io.loadParam(updateGUI=False, rpath=rpath):
            'Error: win::startSimulation load file failed. stop this function.'
            return

    TList=np.linspace(io.T0,io.T1,io.nT)
    HList=np.linspace(io.H0,io.H1,io.nH)
    bondList=[lat.Bond(bond_data[0],bond_data[1],                 # source and target  
                       np.array([int(x) for x in bond_data[2]]),  # over lat.
                       bond_data[3],bond_data[4],bond_data[5],    # strength Jxx, Jyy, Jzz
                       bond_data[6],bond_data[7],bond_data[8],    # strength Jxy, Jxz, Jyz
                       bond_data[9],bond_data[10],bond_data[11],  # strength Jyx, Jzx, Jzy
                       True if io.modelType!='Ising' else False)  # On
                        for bond_data in io.bondList]             # ergodic
    LMatrix=np.array(io.LMatrix)
    pos=np.array(io.pos)
    
    if(io.modelType=='Ising'):
        if io.algorithm!='Metropolis' and io.algorithm!='Wolff':
            print('For now, only Metropolis and Wolff algorithm is supported for Ising model')
            if updateGUI: gui.submitBtn.config(state='normal')
            return
        
        paramPack=[]
        for iH, H in enumerate(HList):
            for iT, T in enumerate(TList):
                paramPack.append([iH*len(TList)+iT,T,bondList,LMatrix,pos,io.S,io.DList,H,io.nsweep,io.nthermal,io.ninterval,io.LPack[0],io.LPack[1],io.LPack[2],io.algorithm,
                                 io.GcOrb,io.orbGroupList,io.groupInSC,io.dipoleAlpha,io.spinFrame])
        
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
        if updateGUI: gui.updateResultViewer(TList=TResult if io.xAxisType=='T' else HResult, magList=[(si+sj)/2 for si,sj in zip(SpinIResult,SpinJResult)], susList=capaResult)
    # continuous model settings
    elif(io.modelType=='XY' or io.modelType=='Heisenberg'):
        for bond in bondList:
            bond.On=True # switch on the vector type bonding
        if io.algorithm!='Metropolis' and io.algorithm!='Wolff':
            print('For now, only Metropolis and Wolff algorithm is supported for O(n) model')
            if updateGUI: gui.submitBtn.config(state='normal')
            return

        On=2 if io.modelType=='XY' else 3
        paramPack=[]
        for iH, H in enumerate(HList):
            for iT, T in enumerate(TList):
                paramPack.append([iH*len(TList)+iT,T,bondList,LMatrix,pos,io.S,io.DList,H,io.nsweep,io.nthermal,io.ninterval,io.LPack[0],io.LPack[1],io.LPack[2],io.algorithm,On,
                                  io.GcOrb,io.orbGroupList,io.groupInSC,io.dipoleAlpha,io.spinFrame])

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
        if updateGUI: gui.updateResultViewer(TList=TResult if io.xAxisType=='T' else HResult, magList=SpinIResult, susList=capaResult)
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