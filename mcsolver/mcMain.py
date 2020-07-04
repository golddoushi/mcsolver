from ctypes import c_double, c_int, CDLL, py_object, c_double
from random import random, randint
import numpy.fft as fft
import Lattice as lat
import numpy as np
import matplotlib.pyplot as plt
import time
import win

class MC:
    def __init__(self,ID,LMatrix,pos=[],S=[],D=[],bondList=[],T=1,Lx=1,Ly=1,Lz=1,ki_s=0,ki_t=0,ki_overLat=[0,0,0],orbGroupList=[],h=0.,dipoleAlpha=0,On=1,spinFrame=0): # init for specified temperature
        self.Lx, self.Ly, self.Lz=Lx, Ly, Lz
        norb=len(pos)
        totOrbs=Lx*Ly*Lz*norb
        # to aviod 0K
        T=0.1 if T<0.1 else T
        # *****************************************************#
        #  ATTENTION: all energies are multiplied by beta here #
        #******************************************************#
        # create orbs for manual temperature
        DT=[d/T for d in D]
        lattice_array, lattice, orbGroup=lat.establishLattice(Lx=Lx,Ly=Ly,Lz=Lz,norb=norb,Lmatrix=np.array(LMatrix),bmatrix=np.array(pos),SpinList=S,DList=DT,orbGroupList=orbGroupList)
        # create bond list for manual temperature
        #bondT=[]
        for bond in bondList:
            bond.renormWithT(T)
            #bondT.append(bond)
        if ki_s>=norb or ki_t>=norb:
            print("ERROR: index out of range ki_S=%d, ki_t=%d, norb=%d\n"%(ki_s,ki_t,norb))
            raise("Input Error!")
            ki_s=norb-1 if ki_s >= norb else ki_s
            ki_t=norb-1 if ki_t >= norb else ki_t
        self.correlatedOrbitalPair=lat.establishLinking(lattice_array,bondList,ki_s=ki_s,ki_t=ki_t,ki_overLat=ki_overLat,dipoleAlpha=dipoleAlpha)
        self.dipoleCorrection = False 
        if abs(dipoleAlpha)>1e-5:
            self.dipoleCorrection=True
            lat.generateDipoleBondings(lattice,dipoleAlpha/T,On=On)
        self.ID=ID
        self.T=T
        self.Sz=Lx*Ly*Lz*sum(S)
        self.Energy=0.
        self.lattice=lattice
        self.lattice_array=lattice_array
        self.orbGroup=orbGroup
        self.totOrbs=totOrbs
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
        initSpin=(c_double*self.totOrbs)()
        nlinking=(c_int*self.totOrbs)()
        nlinking_list=[]
        for iorb, orb in enumerate(self.lattice):
            initSpin[iorb]=c_double(orb.spin)
            nlinking[iorb]=c_int(len(orb.linkedOrb))
            nlinking_list.append(len(orb.linkedOrb))
        
        # link strength
        maxNLinking=np.max(nlinking_list)
        #nlinking=len(orb.linkedOrb)
        linkStrength=(c_double*(self.totOrbs*maxNLinking))()
        #linkStrength_rnorm=(c_double*(self.totOrbs*maxNLinking))() # linking for renormalization
        cnt=0
        for iorb, orb in enumerate(self.lattice):
            for ilinking in range(maxNLinking):
                if ilinking>=nlinking_list[iorb]:
                    linkStrength[cnt]=c_double(0.)
                    #linkStrength_rnorm[cnt]=c_double(0.)
                    cnt+=1
                else:
                    linkStrength[cnt]=c_double(orb.linkStrength[ilinking])
                    #linkStrength_rnorm[cnt]==c_double(orb.linkStrength[ilinking]) if orb.chosen else c_double(0.)
                    cnt+=1
        #for istrength, strength in enumerate(orb.linkStrength):
        #    linkStrength[istrength]=c_double(strength)

        # linking info.
        linkData=(c_int*(self.totOrbs*maxNLinking))()
        cnt=0
        for orb in self.lattice:
            for iorb in range(maxNLinking):#orb.linkedOrb:
                linkData[cnt]=orb.linkedOrb[iorb].id if iorb<len(orb.linkedOrb) else -1
                cnt+=1

        #-------------------------------------------------------------------------------#
        # linking info. for renormalized lattice
        # count total sites in shrinked lat.
        totOrb_rnorm=0
        norbInCluster=-1
        for orb in self.lattice:
            if orb.chosen:
                totOrb_rnorm+=1
                norbInCluster=len(orb.orb_cluster)
        #print("total %d orbs in renormalized lattice"%totOrb_rnorm)
        rOrb=(c_int*totOrb_rnorm)() # store id of renormalized orb cores
        rOrbCluster=(c_int*(totOrb_rnorm*norbInCluster))() # store id of renormalized orbs in cluster
        cnt=0
        for orb in self.lattice:
            if orb.chosen:
                rOrb[cnt]=orb.id
                cnt+=1
        cnt=0
        #print("checking while preparing info.>>>>")
        for orb in self.lattice:
            if orb.chosen:
                #print("orb%d is chosen"%orb.id)
                for orbInCluster in orb.orb_cluster:
                    rOrbCluster[cnt]=orbInCluster.id
                    #print("    orb%d"%orbInCluster.id)
                    cnt+=1
        #print("<<<<")
        
        linkData_rnorm=(c_int*(totOrb_rnorm*maxNLinking))() # store their link info
        cnt=0
        for orb in self.lattice:
            if orb.chosen:
                for iorb in range(maxNLinking):
                    linkData_rnorm[cnt]=c_int(orb.linkedOrb_rnorm[iorb].id) if iorb<len(orb.linkedOrb_rnorm) else c_int(-1)
                    #print("orb%d ~ orb%d, linkeData_rnorm[%d]= %d"%(orb.id,orb.linkedOrb_rnorm[iorb].id,cnt,linkData_rnorm[cnt]))
                    cnt+=1
        #for cnt in range(totOrb_rnorm*maxNLinking):
        #    print(linkData_rnorm[cnt])
        #-------------------------------------------------------------------------------#
        maxNLinking=c_int(maxNLinking)
        
        # field info.
        h=c_double(self.h)
        # correlated info.
        nLat=len(self.correlatedOrbitalPair)
        corrOrbitalPair=(c_int*(nLat*2))()
        for ipair,pair in enumerate(self.correlatedOrbitalPair):
            corrOrbitalPair[ipair*2]=pair[0]
            corrOrbitalPair[ipair*2+1]=pair[1]

        # renormalization calc. switch
        renormOn = c_int(0) if self.dipoleCorrection else c_int(1)

        # num. of returned spin frame
        spinFrame=c_int(int(self.spinFrame))

        mylib=CDLL(win.path+"isinglib.so")
        cMC=mylib.localUpdateMC # default Ising solver
        if algo=='Wolff': # Other choices
            if abs(self.h)>1e-5:
                print("WARNING: external field has not been developed in Wolff algorithm (not effiecient), now the job continues but field is skipped, please switch to Metropolis algorithm to include field effects")
            cMC=mylib.blockUpdateMC
        # uncode returned Data
        cMC.restype=py_object
        data=cMC(self.totOrbs, initSpin, nthermal, nsweep, maxNLinking, nlinking, linkStrength, linkData, ninterval, nLat, corrOrbitalPair, h,
                 renormOn, totOrb_rnorm, norbInCluster, rOrb, rOrbCluster, linkData_rnorm,
                 spinFrame)
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
        self.nsweep=nsweep
        self.nthermal=nthermal
        ninterval=self.totOrbs if ninterval<=0 else ninterval

        # initial spin, single ion anisotropy and number of linking
        initSpin=(c_double*self.totOrbs)()
        initD=(c_double*(3*self.totOrbs))()
        nlinking=(c_int*self.totOrbs)()
        nlinking_list=[]
        for iorb, orb in enumerate(self.lattice):
            #print(orb.spin)
            initSpin[iorb]=c_double(orb.spin)
            initD[iorb*3]=c_double(orb.D[0])
            initD[iorb*3+1]=c_double(orb.D[1])
            initD[iorb*3+2]=c_double(orb.D[2])

            nlinking[iorb]=c_int(len(orb.linkedOrb))
            nlinking_list.append(len(orb.linkedOrb))
        
        # link strength
        maxNLinking=np.max(nlinking_list)
        #print("maxNLinking=%d"%maxNLinking)
        linkStrength=(c_double*(self.totOrbs*maxNLinking*9))() # thus the nlinking of every orbs are the same
        cnt=0
        for iorb, orb in enumerate(self.lattice):
            #print("orb%d"%orb.id)
            for ilinking in range(maxNLinking):
                if ilinking>=nlinking_list[iorb]:
                    for i in range(9):
                        linkStrength[cnt]=c_double(0.)
                        cnt+=1
                else:
                    #print("link %d :"%ilinking,orb.linkStrength[ilinking])
                    for i in range(9):
                        linkStrength[cnt]=c_double(orb.linkStrength[ilinking][i])
                        cnt+=1

        # linking info.
        linkData=(c_int*(self.totOrbs*maxNLinking))()

        #-------------------------------------------------------------------------------#
        # linking info. for renormalized lattice
        # count total sites in shrinked lat.
        totOrb_rnorm=0
        for orb in self.lattice:
            if orb.chosen:
                totOrb_rnorm+=1
                norbInCluster=len(orb.orb_cluster)
        #print("total %d orbs in renormalized lattice"%totOrb_rnorm)
        rOrb=(c_int*totOrb_rnorm)() # store id of renormalized orbs
        rOrbCluster=(c_int*(totOrb_rnorm*norbInCluster))() # store id of renormalized orbs in cluster
        cnt=0
        for orb in self.lattice:
            if orb.chosen:
                rOrb[cnt]=orb.id
                cnt+=1
        
        cnt=0
        #print("checking while preparing info.>>>>")
        for orb in self.lattice:
            if orb.chosen:
                #print("orb%d is chosen"%orb.id)
                for orbInCluster in orb.orb_cluster:
                    rOrbCluster[cnt]=orbInCluster.id
                    #print("    orb%d"%orbInCluster.id)
                    cnt+=1
        #print("<<<<")
        
        linkData_rnorm=(c_int*(totOrb_rnorm*maxNLinking))() # store their link info
        cnt=0
        for orb in self.lattice:
            if orb.chosen:
                for iorb in range(maxNLinking):
                    linkData_rnorm[cnt]=c_int(orb.linkedOrb_rnorm[iorb].id) if iorb<len(orb.linkedOrb_rnorm) else c_int(-1)
                    #print("orb%d ~ orb%d, linkeData_rnorm[%d]= %d"%(orb.id,orb.linkedOrb_rnorm[iorb].id,cnt,linkData_rnorm[cnt]))
                    cnt+=1
        #for cnt in range(totOrb_rnorm*maxNLinking):
        #    print(linkData_rnorm[cnt])
        #-------------------------------------------------------------------------------#


        # field
        h=c_double(self.h)

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

        # correlated info.
        nLat=len(self.correlatedOrbitalPair)
        corrOrbitalPair=(c_int*(nLat*2))()
        for ipair,pair in enumerate(self.correlatedOrbitalPair):
            corrOrbitalPair[ipair*2]=pair[0]
            corrOrbitalPair[ipair*2+1]=pair[1]
            #print(pair[0],pair[1])
        
        # orb group
        nOrbGroup=len(self.orbGroup)
        maxOrbGroupSize=1
        if nOrbGroup>0:
            maxOrbGroupSize=len(self.orbGroup[0])
        if nOrbGroup>1:
            maxOrbGroupSize=np.max([len(subGroup) for subGroup in self.orbGroup])
        orbGroupList=(c_int*(nOrbGroup*maxOrbGroupSize))()
        cnt=0
        for subGroup in self.orbGroup:
            for iorb in range(maxOrbGroupSize):
                if iorb<len(subGroup):
                    orbGroupList[cnt]=subGroup[iorb].id
                else:
                    orbGroupList[cnt]=-1
                cnt+=1
        nOrbGroup=c_int(nOrbGroup)
        maxOrbGroupSize=c_int(maxOrbGroupSize)

        flunc_=c_double(flunc)
        if On==2:
            mylib=CDLL(win.path+"xylib.so")
        elif On==3:
            mylib=CDLL(win.path+"heisenberglib.so")
        else:
            print("Error: undefined O(n)")
            return
        
        spinFrame=c_int(int(self.spinFrame))

        cMC=mylib.localUpdateMC # default On solver
        if algo=='Wolff':
            cMC=mylib.blockUpdateMC
            
        cMC.restype=py_object
        data = cMC(self.totOrbs, initSpin, initD, nthermal, nsweep, maxNLinking_, nlinking, linkStrength, linkData, ninterval, nLat, 
                   corrOrbitalPair, nOrbGroup, maxOrbGroupSize, orbGroupList,
                   flunc_, h,
                   totOrb_rnorm, norbInCluster, rOrb, rOrbCluster, linkData_rnorm,
                   spinFrame)
        spin_i_x, spin_i_y, spin_i_z, spin_j_x, spin_j_y, spin_j_z, spin_ij, autoCorr, E, E2, U4, spin_i_r_x, spin_i_r_y, spin_i_r_z, spin_j_r_x, spin_j_r_y, spin_j_r_z, spin_ij_r, E_r, E2_r, spin_i_tot_z,spin_j_tot_z,spin_tot_z,spin_i_h,spin_j_h,spin_tot_h, spinDistributionList, spinDotSpinBetweenGroups=data
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
        print('%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.3f %.3f %.3f'%(
               self.T,self.h,spin_i_tot_z,spin_j_tot_z,spin_tot_z,spin_i_h,spin_j_h,spin_tot_h,spin_i_len,spin_j_len,spin_ij,np.dot(spin_i,spin_j),E,E2,       U4,spin_ij_r,np.dot(spin_i_r,spin_j_r),       E_r,      E2_r))
        with open('./out','a') as fout:
            fout.write('%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n'%(
                   self.T,self.h,spin_i_tot_z,spin_j_tot_z,spin_tot_z,spin_i_h,spin_j_h,spin_tot_h,spin_i_len,spin_j_len,spin_ij,np.dot(spin_i,spin_j),E,E2,       U4,spin_ij_r,np.dot(spin_i_r,spin_j_r),       E_r,      E2_r))
        if self.spinFrame==nsweep:self.outputSpinWaveSpetra(spinDistributionList)
        
        if self.spinFrame>0:self.outputSpinDistributionForOn(spinDistributionList)
        
        if len(self.orbGroup)>0:self.outputSpinGroup(spinDotSpinBetweenGroups)
        return spin_i, spin_j, spin_ij, autoCorr, E, E2, U4
        
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
            


