import numpy as np
from matplotlib.figure import Figure

class Orbital:
    '''
    represent the orbital, it is linked to many other orbitals
    '''
    def __init__(self,id,spin=1.,D=[],x=0.,y=0.,z=0.,R=np.array([0,0,0])):
        self.id=id
        self.linkedOrb=[]
        self.linkStrength=[]
        self.spin=spin
        self.D=D
        self.inBlock=False

        #use to plot on screen
        self.x=x
        self.y=y
        self.z=z

        self.r=np.array([x,y,z])

        # mark for renormalization
        self.chosen=False
        self.orb_cluster=[]
        self.linkedOrb_rnorm=[]
        self.linkStrength_rnorm=[]
    
    def addLinking(self,targetOrb,strength,quiet=False,forceAdd=False):
        # check redundancy
        for iorb, orb in enumerate(self.linkedOrb):
            if forceAdd:
                break
            if targetOrb.id==orb.id:
                if not quiet: print('Warning: maybe redundant bonding between orb: %d and %d'%(self.id, orb.id))
                # check if they are equal
                existed_bond_strength=self.linkStrength[iorb]
                diff_bond_strength=abs(existed_bond_strength-strength)
                if type(strength) is float and diff_bond_strength<1e-5:
                    if not quiet: print("the difference between this bond and existed bond is negligible (<1e-5), therefore we skip it")
                    return
                if type(strength) is np.ndarray and sum(diff_bond_strength)<1e-5:
                    if not quiet: print("the difference between this bond and existed bond is negligible (<1e-5), therefore we skip it")
                    return
                if not quiet: print('Since the two bond strength is different, now try to add the strength to existed one')
                self.linkStrength[iorb]+=strength
                return
        self.linkedOrb.append(targetOrb)
        self.linkStrength.append(strength)
        #print(strength)
        #exit()
    
    def addLinking_rnorm(self,targetOrb,strength):
        for orb in self.linkedOrb_rnorm:
            if targetOrb.id==orb.id:
                print('Warning: redundant renormalizing bonding between orb: %d and %d'%(self.id, orb.id))
                return
        self.linkedOrb_rnorm.append(targetOrb)
        self.linkStrength_rnorm.append(strength)

    def classifyTheLinking(self,On=False):
        initialType=-1
        self.classStrength=[]
        self.linkedOrbType=[]
        for linkStrength in self.linkStrength:
            #print(linkStrength)
            findType=False
            for itype, StrengthType in enumerate(self.classStrength):
                condition=sum(abs(linkStrength-StrengthType)) if On else abs(linkStrength-StrengthType)
                if condition<0.0001:#abs(linkStrength-StrengthType).all()<0.0001:
                    self.linkedOrbType.append(itype)
                    findType=True
                    break
            if not findType:
                initialType+=1
                self.linkedOrbType.append(initialType)
                self.classStrength.append(linkStrength)
        self.totLinkingTypes=initialType+1

    def getCorrEnergy(self,corrList=[]):
        corr=0.
        for targetOrb, bondStrengh in zip(self.linkedOrb,self.linkStrength):
            excluded=True
            for corrOrb in corrList:
                if targetOrb.id==corrOrb.id:
                    excluded=False
                    break
            if excluded:
                continue
            corr+=self.spin*targetOrb.spin*bondStrengh
        return corr

    def getCorrEnergyDirect(self):
        corr=0.
        for targetOrb, bondStrengh in zip(self.linkedOrb,self.linkStrength):
            corr+=self.spin*targetOrb.spin*bondStrengh
        return corr
    
    def getCorrEnergyWithBlock(self):
        corr=0.
        for targetOrb, bondStrengh in zip(self.linkedOrb,self.linkStrength):
            if targetOrb.inBlock:
                corr+=self.spin*targetOrb.spin*bondStrengh
        return corr

    def addOrbIntoCluster(self,orb_trial):
        newOrb=True
        for orb in self.orb_cluster:
            if orb.id==orb_trial.id:
                newOrb=False
                break
        if newOrb:
            self.orb_cluster.append(orb_trial)
        return
    
class Bond:
    '''
    represent the bond
    '''
    def __init__(self,source,target,overLat,Jxx,Jyy=0,Jzz=0,Jxy=0,Jxz=0,Jyz=0,Jyx=0,Jzx=0,Jzy=0,On=False):
        self.source=source
        self.target=target
        self.overLat=overLat
        self.strength=Jxx
        self.invStrength=Jxx

        self.On=On
        if On:
            self.strength=np.array([Jxx,Jyy,Jzz,Jxy,Jxz,Jyz,Jyx,Jzx,Jzy])
            self.invStrength=np.array([Jxx,Jyy,Jzz,Jyx,Jzx,Jzy,Jxy,Jxz,Jyz])
        #print('new bond')
        #print(self.invStrength)
    
    def renormWithT(self,T):
        if T<1e-4:
            print("Error: Lattice::Bond::renormWithT: Temperature is too low (<1e-4)")
            exit()
        #print(T,self.strength)
        self.strength=(1/T)*self.strength
        self.invStrength=(1/T)*self.invStrength

    def copy(self):
        bond=Bond(self.source,self.target,self.overLat,0,0,0,self.On)
        bond.strength=np.array(list(self.strength)) if self.On else self.strength
        bond.invStrength=np.array(list(self.invStrength)) if self.On else self.invStrength
        return bond

def establishLattice(Lx=1,Ly=1,Lz=1,norb=1,Lmatrix=np.array([[1,0,0],[0,1,0],[0,0,1]]),bmatrix=[np.array([0.,0.,0.])],SpinList=[1],DList=[0.,0.,0.]):
    '''
    create a Lx X Ly X Lz lattice, and create norb orbitals
    for each cell
    '''
    # pre-checking if bmatrix is not consistent with norb
    if len(bmatrix)<norb:
        print('Error: when establish lattice list, we find there is no enough bshift for each orbital')
        exit()

    # now let us begin
    lattice_flatten=[]
    lattice=[]
    id=0
    for x in range(Lx):
        lattice_x=[]
        for y in range(Ly):
            lattice_y=[]
            for z in range(Lz):
                lattice_z=[]
                for o in range(norb):
                    pos=np.dot(np.array([x,y,z])+bmatrix[o],Lmatrix)
                    orbital=Orbital(id,spin=SpinList[o],D=DList[o],
                                    x=pos[0],y=pos[1],z=pos[2],R=np.array([x,y,z]))
                    lattice_z.append(orbital)
                    lattice_flatten.append(orbital)
                    id+=1
                    if x%2+y%2+z%2==0: # mark 1/8 or 1/4 or half orb. for renormalization
                        orbital.chosen=True
                lattice_y.append(lattice_z)
            lattice_x.append(lattice_y)
        lattice.append(lattice_x)
    
    # construct orb cluster for renormalizations
    for x in range(Lx):
        for y in range(Ly):
            for z in range(Lz):
                for o in range(norb):
                    orbital=lattice[x][y][z][o]
                    if orbital.chosen: # 8 possible orbs to add, totally
                        orbital.addOrbIntoCluster(lattice[x][y][z][o])
                        orbital.addOrbIntoCluster(lattice[x][y][(z+1)%Lz][o])
                        orbital.addOrbIntoCluster(lattice[x][(y+1)%Ly][z][o])
                        orbital.addOrbIntoCluster(lattice[(x+1)%Lx][y][z][o])
                        orbital.addOrbIntoCluster(lattice[x][(y+1)%Ly][(z+1)%Lz][o])
                        orbital.addOrbIntoCluster(lattice[(x+1)%Lx][y][(z+1)%Lz][o])
                        orbital.addOrbIntoCluster(lattice[(x+1)%Lx][(y+1)%Ly][z][o])
                        orbital.addOrbIntoCluster(lattice[(x+1)%Lx][(y+1)%Ly][(z+1)%Lz][o])

    # check cluster
    '''print("checking orb cluster after building >>>>>>")
    for orb in lattice_flatten:
        if orb.chosen:
            print("orb%d is chosen, involving:"%orb.id)
            for sub_orb in orb.orb_cluster:
                print("    orb%d"%sub_orb.id)        
    print("<<<<<<")'''       
    return lattice, lattice_flatten

def establishLinking(lattice,bondList,ki_s=0,ki_t=0,ki_overLat=[0,0,0],dipoleAlpha=0):
    Lx=len(lattice)
    Ly=len(lattice[0])
    Lz=len(lattice[0][0])
    Lo=len(lattice[0][0][0])

    correlatedOrbitalPair=[]
    # uncode every orbitals
    for x in range(Lx):
        for y in range(Ly):
            for z in range(Lz):
                for o in range(Lo):
                    # start linking type1: normal bond
                    sourceOrb=lattice[x][y][z][o]
                    for bond in bondList:
                        if o==bond.source:
                            targetOrb=lattice[(x+bond.overLat[0])%Lx][(y+bond.overLat[1])%Ly][(z+bond.overLat[2])%Lz][bond.target]
                            sourceOrb.addLinking(targetOrb,bond.strength)
                            if sourceOrb.id!=targetOrb.id:
                                targetOrb.addLinking(sourceOrb,bond.invStrength)
                    # type2: bond in renormalized system
                    
                    if sourceOrb.chosen:
                        for bond in bondList:
                            if o==bond.source:
                                targetOrb=lattice[(x+bond.overLat[0]*2)%Lx][(y+bond.overLat[1]*2)%Ly][(z+bond.overLat[2]*2)%Lz][bond.target]
                                sourceOrb.addLinking_rnorm(targetOrb,bond.strength)
                                if sourceOrb.id!=targetOrb.id:
                                    targetOrb.addLinking_rnorm(sourceOrb,bond.invStrength)
                # save the correlated orbital pairs
                correlatedOrbitalPair.append([lattice[x][y][z][ki_s].id, lattice[(x+ki_overLat[0])%Lx][(y+ki_overLat[1])%Ly][(z+ki_overLat[2])%Lz][ki_t].id])
    
    # after process
    '''
    On=bondList[0].On
    for x in range(Lx):
        for y in range(Ly):
            for z in range(Lz):
                for o in range(Lo):
                    lattice[x][y][z][o].classifyTheLinking(On=On)
    '''
    return correlatedOrbitalPair

def generateDipoleBondings(lattice,dipoleAlpha,On=1):
    # long-range dipole coupling
    print("Dipole interaction factor %.6f larger than 1e-5, try to construct dipole couplings, note open border condition is employed for dipole interactions"%dipoleAlpha)
    for sourceOrb in lattice:
        for targetOrb in lattice:
            if sourceOrb.id==targetOrb.id:
                continue
            r12=sourceOrb.r-targetOrb.r
            r12_len=np.sqrt(np.dot(r12,r12))
            r12_n=r12/r12_len
            if On==1: # Ising type only AFM part
                dipole_AFM=dipoleAlpha/r12_len**3
                sourceOrb.addLinking(targetOrb,dipole_AFM,forceAdd=True)
                continue
            x,y,z=r12_n
            dipole_AFM=(dipoleAlpha/r12_len**3)*np.eye(3)
            dipole_FM=(-3*dipoleAlpha/r12_len**3)*np.array([[x*x,x*y,x*z],
                                                            [y*x,y*y,y*z],
                                                            [z*x,z*y,z*z]])
            sourceOrb.addLinking(targetOrb,dipole_AFM+dipole_FM,forceAdd=True)
            print("Error reported by Lattice.py::establishLinking  Now dipole interaction for On(2/3) model consists of non-diagonal coupling elements, which is still in development. Therefore, let's stop here.")
            exit()
    print("dipole coupling established")

def plotLattice(lattice):
    '''
    uncode the lattice pack and print each orbital on screen
    '''
    f=Figure()
    ax=f.add_subplot(111)
    for x in lattice:
        for y in x:
            for z in y:
                for o in z:
                    ax.scatter(o.x,o.y,color='blue')
                    ax.annotate(o.id,(o.x,o.y),size=10.5)
                    for target in o.linkedSite:
                        ax.plot([o.x,target.x],[o.y,target.y],color='red')
    f.show()

