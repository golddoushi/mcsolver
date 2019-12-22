from ctypes import c_float, c_int, CDLL, py_object
from random import random, randint
import Lattice as lat
import numpy as np
import time

class MC:
    def __init__(self,ID,LMatrix,pos=[],S=[],D=[],bondList=[],T=1,Lx=1,Ly=1,Lz=1): # init for specified temperature
        norb=len(pos)
        totOrbs=Lx*Ly*Lz*norb
        lattice_array, lattice=lat.establishLattice(Lx=Lx,Ly=Ly,Lz=Lz,norb=norb,Lmatrix=np.array(LMatrix),bmatrix=np.array(pos),SpinList=S,DList=D)
        # to aviod 0K
        T=0.1 if T<0.1 else T
        # create bond list for manual temperature
        bondT=[]
        for bond in bondList:
            bond_tmp=bond.copy()#lat.Bond(bond.source,bond.target,bond.overLat,bond.strength/T)
            bond_tmp.strength/=T;bond_tmp.strength1/=T;bond_tmp.strength2/=T
            bondT.append(bond_tmp)
        lat.establishLinking(lattice_array,bondT)
        self.ID=ID
        self.T=T
        self.Sz=Lx*Ly*Lz*sum(S)
        self.Energy=0.
        self.lattice=lattice
        self.totOrbs=totOrbs
        self.blockLen=0

    def mainLoopViaCLib(self,nsweep=1000,nthermal=5000,algo='Wolff'):
        self.nsweep=nsweep
        self.nthermal=nthermal

        # initial spin
        initSpin=(c_float*self.totOrbs)()
        nlinking=(c_int*self.totOrbs)()
        nlinking_list=[]
        for iorb, orb in enumerate(self.lattice):
            initSpin[iorb]=c_float(orb.spin)
            nlinking[iorb]=c_int(len(orb.linkedOrb))
            nlinking_list.append(len(orb.linkedOrb))
        
        # link strength
        maxNLinking=np.max(nlinking_list)
        #nlinking=len(orb.linkedOrb)
        linkStrength=(c_float*(self.totOrbs*maxNLinking))()
        cnt=0
        for iorb, orb in enumerate(self.lattice):
            for ilinking in range(maxNLinking):
                if ilinking>=nlinking_list[iorb]:
                    linkStrength[cnt]=c_float(0.)
                    cnt+=1
                else:
                    linkStrength[cnt]=c_float(orb.linkStrength[ilinking])
                    cnt+=1
        #for istrength, strength in enumerate(orb.linkStrength):
        #    linkStrength[istrength]=c_float(strength)

        # linking info.
        linkData=(c_int*(self.totOrbs*maxNLinking))()
        cnt=0
        for orb in self.lattice:
            for linkedOrb in orb.linkedOrb:
                linkData[cnt]=linkedOrb.id
                cnt+=1
        
        maxNLinking=c_int(maxNLinking)
        mylib=CDLL("isinglib.so")
        if algo=='Wolff':
            cMC=mylib.blockUpdateMC
            cMC.restype=py_object
            data=cMC(self.totOrbs, initSpin, nthermal, nsweep, maxNLinking, nlinking, linkStrength, linkData)
            print(np.mean(np.abs(data[0]))/self.totOrbs,np.mean(data[1])/self.totOrbs)
            return data[0], data[1]
        elif algo=='Metroplis':
            cMC=mylib.localUpdateMC
            cMC.restype=py_object
            data=cMC(self.totOrbs, initSpin, nthermal, nsweep, maxNLinking, nlinking, linkStrength, linkData)
            print(np.mean(np.abs(data[0]))/self.totOrbs,np.mean(data[1])/self.totOrbs)
            return data[0], data[1]

    def mainLoopViaCLib_On(self,nsweep=1000,nthermal=5000,algo='Metroplis',On=3,flunc=0.0):
        self.nsweep=nsweep
        self.nthermal=nthermal

        # initial spin, single ion anisotropy and number of linking
        initSpin=(c_float*self.totOrbs)()
        initD=(c_float*(3*self.totOrbs))()
        nlinking=(c_int*self.totOrbs)()
        nlinking_list=[]
        for iorb, orb in enumerate(self.lattice):
            #print(orb.spin)
            initSpin[iorb]=c_float(orb.spin)
            initD[iorb*3]=c_float(orb.D[0])
            initD[iorb*3+1]=c_float(orb.D[1])
            initD[iorb*3+2]=c_float(orb.D[2])

            nlinking[iorb]=c_int(len(orb.linkedOrb))
            nlinking_list.append(len(orb.linkedOrb))
        
        # link strength
        maxNLinking=np.max(nlinking_list)
        linkStrength=(c_float*(self.totOrbs*maxNLinking*3))() # thus the nlinking of every orbs are the same
        cnt=0
        for iorb, orb in enumerate(self.lattice):
            for ilinking in range(maxNLinking):
                if ilinking>=nlinking_list[iorb]:
                    for i in range(3):
                        linkStrength[cnt]=c_float(0.)
                        cnt+=1
                else:
                    for i in range(3):
                        linkStrength[cnt]=c_float(orb.linkStrength[ilinking][i])
                        cnt+=1

        # linking info.
        linkData=(c_int*(self.totOrbs*maxNLinking))()
        cnt=0
        for iorb, orb in enumerate(self.lattice):
            for ilinking in range(maxNLinking):
                if ilinking>=nlinking_list[iorb]:
                    linkData[cnt]=-1
                    cnt+=1
                else:
                    linkData[cnt]=orb.linkedOrb[ilinking].id
                    cnt+=1
        
        maxNLinking_=c_int(maxNLinking)
        flunc_=c_float(flunc)
        if On==2:
            mylib=CDLL("xylib.so")
        elif On==3:
            mylib=CDLL("heisenberglib.so")
        else:
            print("Error: undefined O(n)")
            return
        
        if algo=='Wolff':
            #print('Wolff')
            cMC=mylib.blockUpdateMC
            cMC.restype=py_object
            xspin, yspin, zspin, energy = cMC(self.totOrbs, initSpin, initD, nthermal, nsweep, maxNLinking_, nlinking, linkStrength, linkData, flunc_)
            #spin, energy = data[0], data[1], data[2], data[3]
            spin=np.sqrt(np.array(xspin)**2+np.array(yspin)**2+np.array(zspin)**2)
            print('<x> %.3f <y> %.3f <z> %.3f <tot> %.3f <energy> %.3f'%(np.mean(xspin)/self.totOrbs,np.mean(yspin)/self.totOrbs,np.mean(zspin)/self.totOrbs,np.mean(spin)/self.totOrbs,np.mean(energy)/self.totOrbs))
            return spin, np.array(energy)
        elif algo=='Metroplis':
            cMC=mylib.localUpdateMC
            cMC.restype=py_object
            xspin, yspin, zspin, energy = cMC(self.totOrbs, initSpin, initD, nthermal, nsweep, maxNLinking_, nlinking, linkStrength, linkData, flunc_)
            spin = np.sqrt(np.array(xspin)**2+np.array(yspin)**2+np.array(zspin)**2)
            print('<x> %.3f <y> %.3f <z> %.3f <tot> %.3f <energy> %.3f'%(np.mean(np.abs(xspin))/self.totOrbs,np.mean(np.abs(yspin))/self.totOrbs,np.mean(np.abs(zspin))/self.totOrbs,np.mean(np.abs(spin))/self.totOrbs,np.mean(energy)/self.totOrbs))
            return spin, np.array(energy)
        
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
        seedOrb=self.lattice[randint(0,self.totOrbs-1)]
        corr=seedOrb.getCorrEnergyDirect()
        if corr>=0:
            seedOrb.spin*=-1
            self.Sz+=(seedOrb.spin*2)
            self.Energy-=corr*2
        elif np.exp(2*corr)>random():
            seedOrb.spin*=-1
            self.Sz+=(seedOrb.spin*2)
            self.Energy-=corr*2
        return

    def BlockUpdate(self):
        seedOrb=self.lattice[randint(0,self.totOrbs-1)]
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
            


