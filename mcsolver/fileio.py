from tkinter import filedialog
from re import findall
try:
    from . import guiMain as gui
    from . import win
except:
    import guiMain as gui
    import win

global LMatrix, LPack, pos, S, DList, h, H0, H1, nH, dipoleAlpha, bondList, T0, T1, nT, nthermal, nsweep, ninterval, xAxisType, modelType, algorithm, GcOrb, ncores, spinFrame
global orbGroupList, groupInSC
global localCircuitList # local circuit used to calculate the local topological charge
# initial value
xAxisType='T'
orbGroupList=[]
groupInSC=True
GcOrb=[0,0,[0,0,0]]
h=0
H0,H1,nH=0,0,1
dipoleAlpha=0
spinFrame=0
localCircuitList=[]

def collectParam():
    global LMatrix, LPack, pos, S, DList, h, H0, H1, nH, dipoleAlpha, bondList, T0, T1, nT, nthermal, nsweep, ninterval, xAxisType, modelType, algorithm, GcOrb, ncores, spinFrame, localCircuitList
    # get lattice
    a1=gui.latticeGui[0].report()
    a2=gui.latticeGui[1].report()
    a3=gui.latticeGui[2].report()
    LMatrix=[a1,a2,a3]
    print('Lattice matrix:')
    print(a1)
    print(a2)
    print(a3)

    # get supercell size
    LPack=[int(x) for x in gui.supercellGui.report()]
    print('supercell:')
    print(LPack[0],LPack[1],LPack[2])

    # get oribtal position and spin state and onsite-anisotropy
    pos=[ele[3] for ele in gui.OrbListBox.infoData]
    S=[ele[2] for ele in gui.OrbListBox.infoData]
    DList=[ele[4] for ele in gui.OrbListBox.infoData]
    for ipos, iS, iD in zip(pos,S,DList):
        print('positions:',ipos,'Spin:',iS,'onsite-Anisotropy:',iD)

    # set field
    print('isotropic magnetic field is set to %.3f'%h)

    # get bonds
    bondList=[ 
              [bond_data[2][0],bond_data[2][1],\
               bond_data[2][2],\
               bond_data[1][0],bond_data[1][1],bond_data[1][2],bond_data[1][3],bond_data[1][4],bond_data[1][5],bond_data[1][6],bond_data[1][7],bond_data[1][8]]
               for bond_data in gui.BondBox.infoData
             ]
        
    print('bonds:')
    for ibond, bond in enumerate(bondList):
        print("ID %d: orb%d-orb%d [%d %d %d] J:"%(ibond,bond[0],bond[1],bond[2][0],bond[2][1],bond[2][2]))
        print(bond[3:])
    #print(bondList)

    # get TList
    T0, T1, nT=gui.TListGui.report()
    nT=int(nT)
    print('Temperature range: %.2f ~ %.2f with %d sampling points'%(T0, T1, nT))

    # get HList
    H0, H1, nH=gui.HListGui.report()
    nH=int(nH)
    print('Magnetic field range: %.2f ~ %.2f with %d sampling points'%(H0, H1, nH))

    # topological info.
    print('Local circuits per cell: %d'%len(localCircuitList))
    for icircuit, circuit in enumerate(localCircuitList):
        print("ID %d, enclosed by orb %d [%d %d %d], orb %d [%d %d %d], and orb %d [%d %d %d]"%(icircuit,circuit[0][0],*circuit[0][1],circuit[1][0],*circuit[1][1],circuit[2][0],*circuit[2][1]))

    # get thermalizations and sweeps
    nthermal, nsweep, ninterval= [int(x) for x in gui.MCparamGui.report()]
    print('thermalizations, sweeps and tau:')
    print(nthermal, nsweep, ninterval)

    # get model and algorithm
    xAxisType=gui.xAxisGui.get()
    print('X Axis is',xAxisType)
    modelType = gui.modelGui.get()
    print('Model:',modelType)
    algorithm = gui.algorithmGui.get()
    print('Algorithm:',algorithm)

    # get orb. info. for Gc calc.
    #print(gui.corrGui.report())
    s, t, v1, v2, v3 = [int(x) for x in gui.corrGui.report()]
    GcOrb=[[s,t],[v1,v2,v3]]
    print('Measure correlation between orb%d and orb%d with overLat: (%d, %d, %d)'%(s,t,v1,v2,v3))

    # get num. of output spin frames
    spinFrame=int(gui.spinFrameGui.report()[0])
    print('For each setting, %d frames of spin distribution in real space would be output'%spinFrame)

    # get ncores
    ncores= int(gui.coreGui.report()[0])
    print('using %d cores'%ncores)

def saveParam():
    global LMatrix, LPack, pos, S, DList, h, H0, H1, nH, dipoleAlpha, bondList, T0, T1, nT, nthermal, nsweep, ninterval, xAxisType, modelType, algorithm, GcOrb, ncores, spinFrame
    collectParam()
    # write into files
    filePath=filedialog.asksaveasfilename()
    f=open(filePath,'w')
    f.write("This is mcsolver's save file, version: 3.0\n")
    f.write("Lattice:\n")
    a1, a2, a3= LMatrix
    f.write("%.9f %.9f %.9f\n"%(a1[0],a1[1],a1[2]))
    f.write("%.9f %.9f %.9f\n"%(a2[0],a2[1],a2[2]))
    f.write("%.9f %.9f %.9f\n"%(a3[0],a3[1],a3[2]))
    f.write("Supercell used in MC simulations:\n")
    Lx, Ly, Lz=LPack
    f.write("%d %d %d\n"%(Lx,Ly,Lz))
    f.write("Orbitals in cell:\n")
    f.write("%d\n"%len(pos))
    f.write("Positions, initial spin states and onsite-anisotropy of every orbital:\n")
    for ele in gui.OrbListBox.infoData:
        f.write("orb %d: type %d spin %.9f pos [%.9f %.9f %.9f] Dx %.9f Dy %.9f Dz %.9f h %.9f\n"%(ele[0],ele[1],ele[2],\
                                                                                            ele[3][0],ele[3][1],ele[3][2],\
                                                                                            ele[4][0],ele[4][1],ele[4][2],h))
    f.write("Bonds:\n")
    f.write("%d\n"%len(bondList))
    f.write("id, source, target, overLat, exchange matrix elements of each bond:\n")
    for bond_data in gui.BondBox.infoData:
        f.write("bond %d: Jx %.9f Jy %.9f Jz %.9f Jxy %.9f Jxz %.9f Jyz %.9f Jyx %.9f Jzx %.9f Jzy %.9f orb %d to orb %d over [%d %d %d]\n"%\
            (bond_data[0],\
             bond_data[1][0],bond_data[1][1],bond_data[1][2],bond_data[1][3],bond_data[1][4],bond_data[1][5],bond_data[1][6],bond_data[1][7],bond_data[1][8],\
             bond_data[2][0],bond_data[2][1],bond_data[2][2][0],bond_data[2][2][1],bond_data[2][2][2]\
            ))

    f.write("Temperature scanning region:\n")
    f.write("Tmin %.9f Tmax %.9f nT %d\n"%(T0, T1, nT))
    f.write("Field scanning region (in unit 1.48872 T, only if Kelvin and uB is used for energy and spin):\n")
    f.write("Hmin %.9f Hmax %.9f nH %d\n"%(H0, H1, nH))
    f.write("Dipole long-range coupling:\n")
    f.write("alpha %.6f\n"%dipoleAlpha)
    f.write("Measurement:\n")
    f.write("measure the correlation function between orb%d and orb%d over [%d %d %d]\n"%(GcOrb[0][0],GcOrb[0][1],GcOrb[1][0],GcOrb[1][1],GcOrb[1][2]))
    f.write("Supergroup\n")
    f.write("OrbGroup:1\n")
    f.write("Supergroup\n")
    f.write("group0 orb0-orb0\n")
    f.write(">>>       Topological section      <<<\n")
    f.write("LocalCircuit per cell: %d (set to 0 to skip the calc. for topo. Q)\n"%len(localCircuitList))
    for icircuit, circuit in enumerate(localCircuitList):
        f.write("Circuit %d enclosed by orb %d [%d %d %d], orb %d [%d %d %d], and orb %d [%d %d %d]\n"%(icircuit,circuit[0][0],*circuit[0][1],circuit[1][0],*circuit[1][1],circuit[2][0],*circuit[2][1]))
    f.write(">>>   End of Topological section   <<<\n")
    f.write("Distribution output frame: %d\n"%spinFrame)
    f.write("Sweeps for thermalization and statistics, and relaxiation step for each sweep:\n")
    f.write("%d %d %d\n"%(nthermal, nsweep, ninterval))
    f.write("XAxis type:\n")
    f.write(xAxisType+'\n')
    f.write("Model type:\n")
    f.write(modelType+'\n')
    f.write("Algorithm:\n")
    f.write(algorithm+'\n')
    f.write("Ncores:\n")
    f.write("%d\n"%ncores)
    f.close()

def loadParam(updateGUI=True,rpath='./mcInput'):
    global LMatrix, LPack, pos, S, DList, h, bondList, T0, T1, nT, H0, H1, nH, dipoleAlpha, nthermal, nsweep, ninterval, xAxisType, modelType, algorithm, GcOrb, ncores, spinFrame
    global orbGroupList, groupInSC, localCircuitList
    filePath=filedialog.askopenfilename() if updateGUI else rpath
    f=open(filePath,'r')
    data=[line for line in f.read().split('\n') if line]
    f.close()
    # load version info.
    version=findall(r"[0-9\.]+",data[0])[0]
    if version!=str(win.settingFileVersion):
        print("unknown file or version (only support v%s)"%str(win.settingFileVersion))
        return False
    
    # decide position of each tag
    tagLattice=tagSupercell=tagOrbitals=tagBonds=tagTemperature=tagSweeps=tagModel=tagAlgorithm=tagNcores=tagDipole=tagField=tagTopo=0
    for iline, line in enumerate(data):
        keyword=findall(r"[a-zA-Z]+",line)
        if len(keyword)==0:
            continue
        if keyword[0]=='Lattice':
            tagLattice=iline
        if keyword[0]=='Supercell':
            tagSupercell=iline
        if keyword[0]=='Orbitals':
            tagOrbitals=iline
        if keyword[0]=='Bonds':
            tagBonds=iline
        if keyword[0]=='Temperature':
            tagTemperature=iline
        if keyword[0]=='Field':
            tagField=iline
        if keyword[0]=='Dipole':
            tagDipole=iline
        if keyword[0]=='Sweeps':
            tagSweeps=iline
        if keyword[0]=='Model':
            tagModel=iline
        if keyword[0]=='Algorithm':
            tagAlgorithm=iline
        if keyword[0]=='Measurement':
            tagMesurement=iline
        if keyword[0]=='OrbGroup':
            tagOrbGroup=iline
        if keyword[0]=='Ncores':
            tagNcores=iline
        if keyword[0]=='Distribution':
            tagDistribution=iline
        if keyword[0]=='XAxis':
            tagXAxis=iline
        if keyword[0]=='LocalCircuit':
            tagTopo=iline

    if tagLattice*tagSupercell*tagOrbitals*tagBonds*tagTemperature*tagSweeps*tagModel*tagAlgorithm*tagNcores*tagField*tagDipole*tagDistribution*tagXAxis*tagOrbGroup*tagTopo==0:
        print("cannot find some tags")
        return False
    
    # load lattice
    LMatrix=[[float(x) for x in findall(r"[0-9\.\-]+",data[tagLattice+1+i])] for i in range(3)]
    # load supercell
    LPack=[int(x) for x in findall(r"[0-9]+",data[tagSupercell+1])]
    # load orbitals
    norb=int(findall(r"[0-9]+",data[tagOrbitals+1])[0])
    orbInfo=[]
    pos=[]
    S=[]
    DList=[]
    for i in range(norb):
        ele=findall(r"[0-9\.\-]+",data[tagOrbitals+3+i])
        orbInfo.append([int(ele[0]),int(ele[1]),float(ele[2]),\
                       (float(ele[3]),float(ele[4]),float(ele[5])),\
                       (float(ele[6]),float(ele[7]),float(ele[8]))])
        pos.append([float(ele[3]),float(ele[4]),float(ele[5])])
        S.append(float(ele[2]))
        DList.append([float(ele[6]),float(ele[7]),float(ele[8])])
        
    # load bonds
    nbonds=int(findall(r"[0-9]+",data[tagBonds+1])[0])
    bondInfo=[]
    bondList=[]
    for i in range(nbonds):
        ele=findall(r"[0-9\.\-]+",data[tagBonds+3+i])
        bondInfo.append([int(ele[0]),
                         [float(ele[1]),float(ele[2]),float(ele[3]),float(ele[4]),float(ele[5]),float(ele[6]),float(ele[7]),float(ele[8]),float(ele[9])],
                         [int(ele[10]),int(ele[11]),(int(ele[12]),int(ele[13]),int(ele[14]))]
                        ])
                         # source     target       [overlat.                              ] Jz            Jx            Jy            Jxy           Jxz           Jyz           Jyx           Jzy           Jzx
        bondList.append([int(ele[10]),int(ele[11]),[int(ele[12]),int(ele[13]),int(ele[14])],float(ele[1]),float(ele[2]),float(ele[3]),float(ele[4]),float(ele[5]),float(ele[6]),float(ele[7]),float(ele[8]),float(ele[9])])

    # load mesurements
    GcPack=findall(r'[0-9\-]+',data[tagMesurement+1])
    s, t, v1, v2, v3 = [int(x) for x in GcPack]
    GcOrb=[[s,t],[v1,v2,v3]]

    nOrbGroup=int(findall(r"[0-9]+",data[tagOrbGroup])[0])
    groupInSC=True if data[tagOrbGroup+1]=='Supergroup' else False
    orbGroupList=[]
    for i in range(nOrbGroup):
        _, orbID0, orbID1 = findall(r"[0-9]+",data[tagOrbGroup+i+2])
        orbGroupList.append([j for j in range(int(orbID0),int(orbID1)+1)])

    # load topological circuits
    # load bonds
    nCircuits=int(findall(r"[0-9]+",data[tagTopo])[0])
    localCircuitList=[]
    for icircuit in range(nCircuits):
        ele=[int(x) for x in findall(r"[0-9\-]+",data[tagTopo+1+icircuit])]
                             # S1 id  S1 overLat                S2 id  S2 overLat                S3 id  S3 overLat
        circuit_data_packed=[(ele[1],(ele[2],ele[3],ele[4],)),(ele[5],(ele[6],ele[7],ele[8],)),(ele[9],(ele[10],ele[11],ele[12],)),]
        localCircuitList.append(circuit_data_packed)

    # load other parameters
    Tpack=findall(r"[0-9\.]+",data[tagTemperature+1])
    T0, T1, nT = float(Tpack[0]), float(Tpack[1]), int(Tpack[2])
    Hpack=findall(r'[0-9\.\-]+',data[tagField+1])
    H0, H1, nH = float(Hpack[0]), float(Hpack[1]), int(Hpack[2])
    dipoleAlpha=float(findall(r'[0-9\.\-]+',data[tagDipole+1])[0])
    spinFrame=int(findall(r'[0-9]+',data[tagDistribution])[0])
    nTermSweep=findall(r"[0-9\.\-]+",data[tagSweeps+1])
    nthermal, nsweep, ninterval=[int(x) for x in nTermSweep]
    xAxisType=data[tagXAxis+1]
    modelType=data[tagModel+1]
    algorithm=data[tagAlgorithm+1]
    _ncores=[int(data[tagNcores+1])]
    ncores=_ncores[0]

    if not updateGUI:
        return True
    # update gui window
    for i in range(3):
        gui.latticeGui[i].setValue([float(x) for x in LMatrix[i]])
    gui.updateLatticeData()
    gui.supercellGui.setValue(LPack)
    gui.OrbListBox.updateInfo(orbInfo)
    gui.BondBox.updateInfo(bondInfo)
    gui.TListGui.setValue(Tpack)
    gui.HListGui.setValue(Hpack)
    gui.corrGui.setValue(GcPack)
    gui.MCparamGui.setValue(nTermSweep)
    gui.xaxisStr.set(xAxisType)
    gui.modelStr.set(modelType)
    gui.algoStr.set(algorithm)
    gui.spinFrameGui.setValue([spinFrame])
    gui.coreGui.setValue(_ncores)
    gui.updateStructureViewer()
    return True
    