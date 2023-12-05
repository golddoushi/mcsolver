#from ctypes import c_double, c_int, CDLL, py_object, c_double, cdll
from random import random, randint
import numpy.fft as fft
import numpy as np
import matplotlib.pyplot as plt
import time
#import sys,os
try:
    from . import win
    from . import Lattice as lat
except:
    import win
    import Lattice as lat

class MC:
    def __init__(self,ID,LMatrix,pos=[],S=[],D=[],bondList=[],T=1,Lx=1,Ly=1,Lz=1,ki_s=0,ki_t=0,ki_overLat=[0,0,0],orbGroupList=[],groupInSC=False,h=0.,dipoleAlpha=0,On=1,spinFrame=0,localCircuitList=[]): # init for specified temperature
        self.Lx, self.Ly, self.Lz=Lx, Ly, Lz
        norb=len(pos)
        self.totOrbs=Lx*Ly*Lz*norb
        # to aviod 0K
        T=0.1 if T<0.1 else T
        # *****************************************************#
        #  ATTENTION: all energies are multiplied by beta here #
        #******************************************************#
        # create orbs for manual temperature
        DT=[d/T for d in D]
        self.lattice_array, self.lattice, self.orbGroup, self.localCircuit=lat.establishLattice(Lx=Lx,Ly=Ly,Lz=Lz,norb=norb,Lmatrix=np.array(LMatrix),bmatrix=np.array(pos),SpinList=S,DList=DT,orbGroupList=orbGroupList,groupInSC=groupInSC,localCircuitList=localCircuitList)
        # create bond list for manual temperature
        #bondT=[]
        for bond in bondList:
            bond.renormWithT(T)
            #bondT.append(bond)
        if ki_s>=norb or ki_t>=norb:
            print("ERROR: index out of range ki_S=%d, ki_t=%d, norb=%d\n"%(ki_s,ki_t,norb))
            raise("Input Error!")
        self.correlatedOrbitalPair=lat.establishLinking(self.lattice_array,bondList,ki_s=ki_s,ki_t=ki_t,ki_overLat=ki_overLat,Lmatrix=np.array(LMatrix),bmatrix=np.array(pos),dipoleAlpha=dipoleAlpha)
        lat.constructLocalFrame(self.lattice)
        self.dipoleCorrection = False 
        if abs(dipoleAlpha)>1e-5:
            self.dipoleCorrection=True
            lat.generateDipoleBondings(self.lattice,dipoleAlpha/T,On=On)
        self.ID=ID
        self.T=T
        self.Sz=Lx*Ly*Lz*sum(S)
        self.Energy=0.
        self.blockLen=0
        self.ki_s=ki_s
        self.ki_t=ki_t
        self.h=h
        self.spinFrame=spinFrame
    
    def mainLoopViaCLib(self,nsweep=1000,nthermal=500,ninterval=-1,algo='Wolff'):
        self.nsweep=nsweep
        self.nthermal=nthermal
        ninterval=self.totOrbs if ninterval<=0 else ninterval

        # initial spin
        initSpin=[]
        nlinking=[]
        nlinking_list=[]
        for iorb, orb in enumerate(self.lattice):
            initSpin.append(orb.spin)
            nlinking.append(len(orb.linkedOrb))
            nlinking_list.append(len(orb.linkedOrb))
        
        # link strength
        maxNLinking=max(nlinking_list)
        #nlinking=len(orb.linkedOrb)
        linkStrength,linkData=[],[]
        #linkStrength_rnorm=(c_double*(self.totOrbs*maxNLinking))() # linking for renormalization
        for iorb, orb in enumerate(self.lattice):
            for ilinking in range(maxNLinking):
                if ilinking>=nlinking_list[iorb]:
                    linkData.append(-1)
                    linkStrength.append(0.)
                    #linkStrength_rnorm[cnt]=c_double(0.)
                else:
                    linkData.append(orb.linkedOrb[ilinking].id)
                    linkStrength.append(orb.linkStrength[ilinking])

        #-------------------------------------------------------------------------------#
        # linking info. for renormalized lattice
        # count total sites in shrinked lat.
        
        #print("total %d orbs in renormalized lattice"%totOrb_rnorm)
        rOrb,rOrbCluster,linkData_rnorm=[],[],[] # store id of renormalized orbs in cluster
        for orb in self.lattice:
            if orb.chosen:
                rOrb.append(orb.id)
                #print("orb%d is chosen"%orb.id)
                for orbInCluster in orb.orb_cluster:
                    rOrbCluster.append(orbInCluster.id)

                for iorb in range(maxNLinking):
                    if iorb<len(orb.linkedOrb_rnorm):
                        linkData_rnorm.append(orb.linkedOrb_rnorm[iorb].id) 
                    else:
                        linkData_rnorm.append(-1)
                    #print("orb%d ~ orb%d, linkeData_rnorm[%d]= %d"%(orb.id,orb.linkedOrb_rnorm[iorb].id,cnt,linkData_rnorm[cnt]))
        #for cnt in range(totOrb_rnorm*maxNLinking):
        #    print(linkData_rnorm[cnt])
        #-------------------------------------------------------------------------------#
        
        # correlated info.
        corrOrbitalPair=[]
        for pair in self.correlatedOrbitalPair: corrOrbitalPair+=pair

        updateAlgorithm=0 # default Metropolis algorithm
        if algo=='Wolff': updateAlgorithm=1

        try:
            from isinglib import MCMainFunction
        except:
            from mcsolver.lib.isinglib import MCMainFunction

        def callback(k):
            print("callbak function of python is reporting, recived parameter is:")
            print(k)
            return
        data=MCMainFunction(updateAlgorithm, tuple(initSpin), nthermal, nsweep, ninterval, 
                            maxNLinking, tuple(nlinking), tuple(linkStrength), tuple(linkData), 
                            tuple(corrOrbitalPair), self.h/self.T,
                            tuple(rOrb), tuple(rOrbCluster), tuple(linkData_rnorm),
                            self.spinFrame,
                            callback)
        #print('data has been returned successfully, dim=%d'%len(data))
        spin_i, spin_j, spin_ij, autoCorr, E, E2, E_rnorm, E2_rnorm, U4, spin_tot, spinDistributionList = data
        E*=self.T;E2*=self.T**2;E_rnorm*=self.T;E2_rnorm*=self.T**2 # recover the real energies
        print("%.3f %.3f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.6f %.6f %.6f"%(
               self.T, self.h,   spin_i,    spin_j,  spin_tot,   spin_ij,        autoCorr,    E,        E2,   E_rnorm,    E2_rnorm, E2-E**2, E2_rnorm-E_rnorm**2, U4))
        with open('./out','a') as fout:
            fout.write("%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f\n"%(
                self.T, self.h,   spin_i,    spin_j,  spin_tot,   spin_ij,        E,        E2,   E_rnorm,    E2_rnorm, U4))
        if self.spinFrame>0:
            self.outputSpinDistributionForIsing(spinDistributionList)
        return spin_i, spin_j, spin_ij, autoCorr, E, E2, U4

    def mainLoopViaCLib_On(self,nsweep=1000,nthermal=5000,ninterval=-1,algo='Metroplis',On=3,flunc=0.0,h=0.,binGraph=False):
        def callback(k):
            print("callbak function of python is reporting, recived parameter is:")
            print(k)
            return
        
        self.nsweep=nsweep
        self.nthermal=nthermal
        ninterval=self.totOrbs if ninterval<=0 else ninterval

        updateAlgorithm=0 # default Metropolis algorithm
        if algo=='Wolff': updateAlgorithm=1

        initSpin,initD,nlinking=[],[],[]
        for orb in self.lattice:
            initSpin.append(orb.spin)
            for i in range(3):
                initD.append(orb.D[i])
            nlinking.append(len(orb.linkedOrb))
        
        # link strength
        ignoreNonDiagonalJ=1
        maxNLinking=max(nlinking)
        linkStrength,linkData=[],[] # thus the nlinking of every orbs are the same
        
        for iorb, orb in enumerate(self.lattice):
            #print("orb%d"%orb.id)
            for ilinking in range(maxNLinking):
                if ilinking>=nlinking[iorb]:
                    linkData.append(-1)
                    for i in range(9):  # set the redundant bond strength to zero
                        linkStrength.append(0.)
                else:
                    linkData.append(orb.linkedOrb[ilinking].id)
                    #print("link %d :"%ilinking,orb.linkStrength[ilinking])
                    for i in range(9):  # set the bond strength
                        linkStrength.append(orb.linkStrength[ilinking][i])
                        if abs(orb.linkStrength[ilinking][i]) > 1e-6 and i>=3: ignoreNonDiagonalJ=0

        # topological circuits
        localCircuits=[]
        for circuit in self.localCircuit:
            for orb in circuit: localCircuits+=[orb.id]
        
        # correlated info.
        corrOrbitalPair=[]
        for pair in self.correlatedOrbitalPair:corrOrbitalPair+=pair

        # orb group
        nOrbGroup=len(self.orbGroup)
        maxOrbGroupSize=1
        if nOrbGroup>0:
            maxOrbGroupSize=len(self.orbGroup[0])
        if nOrbGroup>1:
            maxOrbGroupSize=np.max([len(subGroup) for subGroup in self.orbGroup])
        orbGroupList=[]
        for subGroup in self.orbGroup:
            for iorb in range(maxOrbGroupSize):
                if iorb<len(subGroup):
                    orbGroupList.append(subGroup[iorb].id)
                else:
                    orbGroupList.append(-1)
        #print(nOrbGroup,maxOrbGroupSize)
        #exit()

        #-------------------------------------------------------------------------------#
        # linking info. for renormalized lattice
        # count total sites in shrinked lat.
        
        #print("total %d orbs in renormalized lattice"%totOrb_rnorm)
        rOrb,rOrbCluster,linkData_rnorm=[],[],[] # store id of renormalized orbs
        #print("checking while preparing info.>>>>")
        for orb in self.lattice:
            if orb.chosen:
                #print("orb%d is chosen"%orb.id)
                rOrb.append(orb.id)
                for orbInCluster in orb.orb_cluster:
                    rOrbCluster.append(orbInCluster.id)
                    #print("    orb%d"%orbInCluster.id)
            
                for iorb in range(maxNLinking):
                    if iorb<len(orb.linkedOrb_rnorm):
                        linkData_rnorm.append(orb.linkedOrb_rnorm[iorb].id)  
                    else:
                        linkData_rnorm.append(-1)
                    #print("orb%d ~ orb%d, linkeData_rnorm[%d]= %d"%(orb.id,orb.linkedOrb_rnorm[iorb].id,cnt,linkData_rnorm[cnt]))
        #print("<<<<")
                
        #for cnt in range(totOrb_rnorm*maxNLinking):
        #    print(linkData_rnorm[cnt])
        #-------------------------------------------------------------------------------#
        if On==2:
            try:
                from xylib import MCMainFunction
            except:
                from mcsolver.lib.xylib import MCMainFunction
        elif On==3:
            try:
                from heisenberglib import MCMainFunction
            except:
                from mcsolver.lib.heisenberglib import MCMainFunction
        
        data = MCMainFunction(updateAlgorithm,tuple(initSpin),tuple(initD),nthermal,nsweep,ninterval,
                             maxNLinking,tuple(nlinking),tuple(linkStrength),tuple(linkData),
                             tuple(localCircuits),
                             tuple(corrOrbitalPair),
                             nOrbGroup, maxOrbGroupSize, tuple(orbGroupList),
                             flunc, self.h/self.T,
                             tuple(rOrb), tuple(rOrbCluster), tuple(linkData_rnorm),
                             self.spinFrame,
                             ignoreNonDiagonalJ,
                             callback)
        
        spin_i_x, spin_i_y, spin_i_z, spin_j_x, spin_j_y, spin_j_z, spin_ij, autoCorr, E, E2, U4, spin_i_r_x, spin_i_r_y, spin_i_r_z, spin_j_r_x, spin_j_r_y, spin_j_r_z, spin_ij_r, E_r, E2_r, spin_i_tot_z,spin_j_tot_z,spin_tot_z,spin_i_h,spin_j_h,spin_tot_h, topologicalQ, spinDistributionList, spinDotSpinBetweenGroups=data
        E*=self.T;E2*=self.T**2;E_r*=self.T;E2_r*=self.T**2 # recover the real energies
        C=E2-E*E
        C_r=E2_r-E_r*E_r
        spin_i=np.array([spin_i_x, spin_i_y, spin_i_z])
        spin_j=np.array([spin_j_x, spin_j_y, spin_j_z])
        spin_i_len=np.sqrt(np.dot(spin_i,spin_i))
        spin_j_len=np.sqrt(np.dot(spin_j,spin_j))
        spin_i_r=np.array([spin_i_r_x, spin_i_r_y, spin_i_r_z])
        spin_j_r=np.array([spin_j_r_x, spin_j_r_y, spin_j_r_z])
        #      T       <i><j>     <ij>      <autoCorr>      <E>      <E2>      <U4>      <E_r>      <E2_r>  C  C_v
        print('T=%.3f h=%.3f <Siz>=%.3f <Sjz>=%.3f <Sz>=%.3f <Sih>=%.3f <Sjh>=%.3f <Sh>=%.3f <Si>=%.3f <Sj>=%.3f <SiSj>=%.3f <Si><Sj>=%.6f <E>=%.3f <E^2>=%.3f <U4>=%.3f <SiSj>r=%.3f <Si>r<Sj>r=%.3f <Er>=%.3f <E^2_r>=%.3f <Q>=%.3f'%(
               self.T,self.h,spin_i_tot_z,spin_j_tot_z,spin_tot_z,spin_i_h,spin_j_h,spin_tot_h,spin_i_len,spin_j_len,spin_ij,np.dot(spin_i,spin_j),E,E2,       U4,spin_ij_r,np.dot(spin_i_r,spin_j_r),       E_r,      E2_r, topologicalQ))
        with open('./out','a') as fout:
            fout.write('T= %.6E h= %.6E <Siz>= %.6E <Sjz>= %.6E <Sz>= %.6E <Sih>= %.6E <Sjh>= %.6E <Sh>= %.6E <Si>= %.6E <Sj>= %.6E <SiSj>= %.6E <Si><Sj>= %.6E <E>= %.6E <E^2>= %.6E <U4>= %.6E <SiSj>r= %.6E <Si>r<Sj>r= %.6E <Er>= %.6E <E^2_r>= %.6E <Q>= %.6E\n'%(
                   self.T,self.h,spin_i_tot_z,spin_j_tot_z,spin_tot_z,spin_i_h,spin_j_h,spin_tot_h,spin_i_len,spin_j_len,spin_ij,np.dot(spin_i,spin_j),E,E2,       U4,spin_ij_r,np.dot(spin_i_r,spin_j_r),       E_r,      E2_r, topologicalQ))
        # FFT on the MC random data makes no sense
        #if self.spinFrame==nsweep:self.outputSpinWaveSpetra(spinDistributionList) 
        
        if self.spinFrame>0:self.outputSpinDistributionForOn(spinDistributionList)
        
        if len(self.orbGroup)>0:self.outputSpinGroup(spinDotSpinBetweenGroups)
        return spin_i, spin_j, spin_ij, autoCorr, E, E2, U4, topologicalQ

    def outputSpinGroup(self,spinDotSpinData):
        with open('./spinDotSpin.txt','a') as fout:
            # title
            fout.write('%.3f %.3f '%(self.T,self.h))
            cnt=0
            for i in range(len(self.orbGroup)+1):
                for j in range(len(self.orbGroup)+1):
                    #keyword1='group_%d'%i if i!=self.orbGroup else 'Total'
                    #keyword2='group_%d'%j if j!=self.orbGroup else 'Total'
                    #title='#'+keyword1+'_'+keyword2+' '
                    fout.write('%.6f '%spinDotSpinData[cnt])
                    cnt+=1
            for i in range(len(self.orbGroup)+1):
                fout.write('%.6f '%spinDotSpinData[cnt])
                cnt+=1
            fout.write('\n')

    def outputSpinDistributionForIsing(self,distributionList):
        for iframe in range(self.spinFrame):
            with open('./IsingSpinDistribution.T%.3f.H%.3f.%d.txt'%(self.T,self.h,iframe),'w') as fout:
                fout.write("#x       #y       #z       #spin\n")
                for orb in self.lattice:
                    fout.write('%.6f %.6f %.6f %.3f\n'%(orb.x,orb.y,orb.z,
                                distributionList[iframe][orb.id]))
    
    def outputSpinWaveSpetra(self,distributionList):
        print("start FFT")
        time0=time.time()
        # get the 1st orb on gamma-K
        cnt=0
        id_list=[]
        for i in range(self.Lx):
            id_list.append(self.lattice_array[i][i][0][0].id)
            cnt+=1
            if cnt>=self.Ly:
                break
        T_R_x_data=np.array(distributionList)[:,id_list,0]
        amplitude=np.log(abs(fft.fft2(T_R_x_data)))[:50,:]
        print("end FFT, time elapsed %.3fs"%(time.time()-time0))
        plt.imshow(amplitude,origin='lower',extent=(0,1,0,1))
        plt.show()
        exit()

    def outputSpinDistributionForOn(self,distributionList):
        for iframe in range(self.spinFrame):
            with open('./OnSpinDistribution.T%.3f.H%.3f.%d.txt'%(self.T,self.h,iframe),'w') as fout:
                fout.write("#x       #y       #z       #spinx  #spiny  #spinz\n")
                for orb in self.lattice:
                    fout.write('%.6f %.6f %.6f %.6f %.6f %.6f\n'%(orb.x,orb.y,orb.z,
                                distributionList[iframe][orb.id][0],distributionList[iframe][orb.id][1],distributionList[iframe][orb.id][2]))

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
            


