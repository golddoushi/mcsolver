import numpy as np
from matplotlib.figure import Figure

class Orbital:
    '''
    represent the orbital, it is linked to many other orbitals
    '''
    def __init__(self,id,spin=1.,x=0.,y=0.,z=0.):
        self.id=id
        self.linkedOrb=[]
        self.linkStrength=[]
        self.spin=spin
        self.inBlock=False

        #use to plot on screen
        self.x=x
        self.y=y
        self.z=z
    
    def addLinking(self,targetOrb,strength):
        #print(strength)
        # check redundancy
        for orb in self.linkedOrb:
            if targetOrb.id==orb.id:
                print('Warning: redundant bonding between orb: %d and %d'%(self.id, orb.id))
                return
        self.linkedOrb.append(targetOrb)
        self.linkStrength.append(strength)

    def classifyTheLinking(self):
        initialType=-1
        self.classStrength=[]
        self.linkedOrbType=[]
        for linkStrength in self.linkStrength:
            findType=False
            for itype, StrengthType in enumerate(self.classStrength):
                if abs(linkStrength-StrengthType).all()<0.0001:
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
    

class Bond:
    '''
    represent the bond
    '''
    def __init__(self,source,target,overLat,strength,strength1=0.,strength2=0.,On=False):
        self.source=source
        self.target=target
        self.overLat=overLat
        self.strength=strength

        self.On=On
        self.strength1=strength1 # for xy and heisenberg model
        self.strength2=strength2 # for xy and heisenberg model
    
    def copy(self):
        return Bond(self.source,self.target,self.overLat,self.strength,self.strength1,self.strength2,self.On)

def establishLattice(Lx=1,Ly=1,Lz=1,norb=1,Lmatrix=np.array([[1,0,0],[0,1,0],[0,0,1]]),bmatrix=[np.array([0.,0.,0.])],SpinList=[1]):
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
                    orbital=Orbital(id,spin=SpinList[o],x=pos[0],y=pos[1],z=pos[2])
                    lattice_z.append(orbital)
                    lattice_flatten.append(orbital)
                    id+=1
                lattice_y.append(lattice_z)
            lattice_x.append(lattice_y)
        lattice.append(lattice_x)
    return lattice, lattice_flatten

def establishLinking(lattice,bondList):
    Lx=len(lattice)
    Ly=len(lattice[0])
    Lz=len(lattice[0][0])
    Lo=len(lattice[0][0][0])
    # uncode every orbitals
    for x in range(Lx):
        for y in range(Ly):
            for z in range(Lz):
                for o in range(Lo):
                    # start linking
                    sourceOrb=lattice[x][y][z][o]
                    for bond in bondList:
                        if o==bond.source:
                            targetOrb=lattice[(x+bond.overLat[0])%Lx][(y+bond.overLat[1])%Ly][(z+bond.overLat[2])%Lz][bond.target]
                            if bond.On:
                                sourceOrb.addLinking(targetOrb,np.array([bond.strength,bond.strength1,bond.strength2]))
                            else:
                                sourceOrb.addLinking(targetOrb,bond.strength)
                            if sourceOrb.id!=targetOrb.id:
                                if bond.On:
                                    targetOrb.addLinking(sourceOrb,np.array([bond.strength,bond.strength1,bond.strength2]))
                                else:
                                    targetOrb.addLinking(sourceOrb,bond.strength)
    # after process
    for x in range(Lx):
        for y in range(Ly):
            for z in range(Lz):
                for o in range(Lo):
                    lattice[x][y][z][o].classifyTheLinking()

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

