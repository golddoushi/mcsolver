import Lattice as lat
import numpy as np
from ctypes import *
import random
import time

class MC:
    def __init__(self,ID,LMatrix,pos,bondList,T=1,Lx=1,Ly=1,Lz=1): # init for specified temperature
        norb=len(pos)
        totOrbs=Lx*Ly*Lz*norb
        lattice_array, lattice=lat.establishLattice(Lx=Lx,Ly=Ly,Lz=Lz,norb=norb,Lmatrix=np.array(LMatrix),bmatrix=np.array(pos))
        # create bond list for manual temperature
        bondT=[]
        for bond in bondList:
            bond_tmp=lat.Bond(bond.source,bond.target,bond.overLat,bond.strength/T)
            bondT.append(bond_tmp)
        lat.establishLinking(lattice_array,bondT)
        self.ID=ID
        self.T=T
        self.Sz=totOrbs
        self.Energy=0.
        self.lattice=lattice
        self.totOrbs=totOrbs
        self.blockLen=0

    def mainLoopViaCLib(self,nsweep=1000,nthermal=5000,algo='Wolff'):
        self.nsweep=nsweep
        self.nthermal=nthermal

        # initial spin
        #SpinField=c_float*self.totOrbs
        initSpin=(c_float*self.totOrbs)()
        for iorb, orb in enumerate(self.lattice):
            initSpin[iorb]=c_float(orb.spin)
        
        # link strength
        nlinking=len(orb.linkedOrb)
        linkStrength=(c_float*nlinking)()
        for istrength, strength in enumerate(orb.linkStrength):
            linkStrength[istrength]=c_float(strength)

        # linking info.
        linkData=(c_int*(self.totOrbs*nlinking))()
        cnt=0
        for orb in self.lattice:
            for linkedOrb in orb.linkedOrb:
                linkData[cnt]=linkedOrb.id
                cnt+=1

        time0=time.time()
        mylib=CDLL("isinglib.so")
        if algo=='Wolff':
            cMC=mylib.blockUpdateMC
            cMC.restype=py_object
            data=cMC(self.totOrbs, initSpin, nthermal, nsweep, nlinking, linkStrength, linkData)
            print(np.mean(np.abs(data[0]))/self.totOrbs,np.mean(data[1])/self.totOrbs,time.time()-time0)
            return data[0], data[1]
        elif algo=='Metroplis':
            cMC=mylib.localUpdateMC
            cMC.restype=py_object
            data=cMC(self.totOrbs, initSpin, nthermal, nsweep, nlinking, linkStrength, linkData)
            print(np.mean(np.abs(data[0]))/self.totOrbs,np.mean(data[1])/self.totOrbs,time.time()-time0)
            return data[0], data[1]
        

    def mainLoop(self,nsweep=10000,nthermal=5000):
        self.nsweep=nsweep
        self.nthermal=nthermal
        sAvgData=[]
        EnergyData=[]
        blockData=[]
    
        process=0
        for ithermal in range(nthermal):
            if ithermal>=nthermal*0.01*process:
                #print('thermalization %2d percent'%process)
                process+=1
            #print('thermalization:',ithermal)
            for imcStep in range(self.totOrbs):
                #self.LocalUpdate()
                self.BlockUpdate()

        print('sweep started')
        process=0
        for isweep in range(nsweep):
            if isweep>=nsweep*0.01*process:
                print('simulation has done %2d percent'%process)
                process+=1
            sAvgSweep=energySweep=blockSweep=0
            for imcStep in range(self.totOrbs):
                #self.LocalUpdate()
                self.BlockUpdate()
                sAvgSweep+=np.abs(self.Sz)
                energySweep+=self.Energy
                blockSweep+=self.blockLen
            
            #print(sAvgSweep/self.totOrbs)
            sAvgData.append(sAvgSweep/self.totOrbs)
            EnergyData.append(energySweep/self.totOrbs)
            blockData.append(blockSweep/self.totOrbs)
        
        self.sAvgData=sAvgData
        self.EnergyData=EnergyData
        print(np.mean(sAvgData)/self.totOrbs, np.std(sAvgData)/self.totOrbs, np.mean(blockData)/self.totOrbs/self.T)
        return np.mean(sAvgData), np.mean(EnergyData), np.std(sAvgData), np.std(EnergyData)

    def saveData(self):
        data=open('./task'+str(self.ID),'w')
        data.write('Temperature:%.3f \n'%self.T)
        data.write('Energy List, %5d total lines:\n'%self.nsweep)
        for energy in self.EnergyData:
            data.write('%.6f \n'%(energy/self.totOrbs))
        data.write('spin List, %5d total lines: \n'%self.nsweep)
        for spin in self.sAvgData:
            data.write('%.6f \n'%abs(spin/self.totOrbs))
        data.close()

    def LocalUpdate(self):
        seedOrb=self.lattice[random.randint(0,self.totOrbs-1)]
        corr=seedOrb.getCorrEnergyDirect()
        if corr>=0:
            seedOrb.spin*=-1
            self.Sz+=(seedOrb.spin*2)
            self.Energy-=corr*2
        elif np.exp(2*corr)>random.random():
            seedOrb.spin*=-1
            self.Sz+=(seedOrb.spin*2)
            self.Energy-=corr*2
        return

    def BlockUpdate(self):
        seedOrb=self.lattice[random.randint(0,self.totOrbs-1)]
        seedOrb.inBlock=True
        block=[seedOrb]
        buffer=[seedOrb]

        def expandBuffer():
            #print(len(buffer))
            if len(buffer)==0:
                return False
            outOrb=buffer.pop(0)
            for ilinkedOrb, linkedOrb in enumerate(outOrb.linkedOrb):
                if linkedOrb.inBlock:
                    continue
                corr=outOrb.linkStrength[ilinkedOrb]*outOrb.spin*linkedOrb.spin
                if corr<0 and (1-np.exp(2*corr))>random.random():
                    linkedOrb.inBlock=True
                    block.append(linkedOrb)
                    buffer.append(linkedOrb)
            return True

        while(expandBuffer()):
            pass

        for blockOrb in block:
            blockOrb.spin*=-1
            blockOrb.inBlock=False
            self.Sz+=(blockOrb.spin*2)
        #print(len(block))
        self.blockLen=len(block)
        return
            


