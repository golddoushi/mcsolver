from tkinter import filedialog
from re import findall
import guiMain as gui

global LMatrix, LPack, pos, S, DList, bondList, T0, T1, nT, nthermal, nsweep, ninterval, modelType, algorithm, GcOrb, ncores
# initial value
GcOrb=[0,0,[0,0,0]]

def collectParam():
    global LMatrix, LPack, pos, S, DList, bondList, T0, T1, nT, nthermal, nsweep, ninterval, modelType, algorithm, GcOrb, ncores
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

    # get bonds
    bondList=[
              [bond_data[2][0],bond_data[2][1],\
               bond_data[2][2],\
               bond_data[1][0],bond_data[1][1],bond_data[1][2]] \
               for bond_data in gui.BondBox.infoData
             ]
        
    print('bonds:')
    print(bondList)

    # get TList
    T0, T1, nT=gui.TListGui.report()
    nT=int(nT)
    print('Temperature range: %.2f ~ %.2f with %d sampling points'%(T0, T1, nT))

    # get thermalizations and sweeps
    nthermal, nsweep, ninterval= [int(x) for x in gui.MCparamGui.report()]
    print('thermalizations, sweeps and tau:')
    print(nthermal, nsweep, ninterval)

    # get model and algorithm
    modelType = gui.modelGui.get()
    print('Model:',modelType)
    algorithm = gui.algorithmGui.get()
    print('Algorithm:',algorithm)

    # get orb. info. for Gc calc.
    #print(gui.corrGui.report())
    s, t, v1, v2, v3 = [int(x) for x in gui.corrGui.report()]
    GcOrb=[[s,t],[v1,v2,v3]]
    print('Measure correlation between orb%d and orb%d with overLat: (%d, %d, %d)'%(s,t,v1,v2,v3))

    # get ncores
    ncores= int(gui.coreGui.report()[0])
    print('using %d cores'%ncores)

def saveParam():
    global LMatrix, LPack, pos, S, DList, bondList, T0, T1, nT, nthermal, nsweep, ninterval, modelType, algorithm, GcOrb, ncores
    collectParam()
    # write into files
    filePath=filedialog.asksaveasfilename()
    f=open(filePath,'w')
    f.write("This is mcsolver's save file, version: 1.2\n")
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
    f.write("Postions, initial spin states and onsite-anisotropy of every orbitals:\n")
    for ele in gui.OrbListBox.infoData:
        f.write("orb %d: type %d spin %.9f pos [%.9f %.9f %.9f] Dz %.9f Dx %.9f Dy %.9f\n"%(ele[0],ele[1],ele[2],\
                                                                                            ele[3][0],ele[3][1],ele[3][2],\
                                                                                            ele[4][0],ele[4][1],ele[4][2]))
    f.write("Bonds:\n")
    f.write("%d\n"%len(bondList))
    f.write("id, source, target, overLat, Jz, Jx, Jy of each bond:\n")
    for bond_data in gui.BondBox.infoData:
        f.write("bond %d: Jz %.9f Jx %.9f Jy %.9f orb %d to orb %d over [%d %d %d]\n"%\
            (bond_data[0],\
             bond_data[1][0],bond_data[1][1],bond_data[1][2],\
             bond_data[2][0],bond_data[2][1],bond_data[2][2][0],bond_data[2][2][1],bond_data[2][2][2]\
            ))

    f.write("Temperature scanning region:\n")
    f.write("Tmin %.9f Tmax %.9f nT %d\n"%(T0, T1, nT))
    f.write("Mesurement:\n")
    f.write("mesure the correlation function between orb%d and orb%d over [%d %d %d]\n"%(GcOrb[0][0],GcOrb[0][1],GcOrb[1][0],GcOrb[1][1],GcOrb[1][2]))
    f.write("Sweeps for thermalization and statistics, and relaxiation step for each sweep:\n")
    f.write("%d %d %d\n"%(nthermal, nsweep, ninterval))
    f.write("Model type:\n")
    f.write(modelType+'\n')
    f.write("Algorithm:\n")
    f.write(algorithm+'\n')
    f.write("Ncores:\n")
    f.write("%d\n"%ncores)
    f.close()

def loadParam(updateGUI=True,rpath='./mcInput'):
    global LMatrix, LPack, pos, S, DList, bondList, T0, T1, nT, nthermal, nsweep, ninterval, modelType, algorithm, GcOrb, ncores
    filePath=filedialog.askopenfilename() if updateGUI else rpath
    f=open(filePath,'r')
    data=[line for line in f.read().split('\n') if line]
    f.close()
    # load version info.
    version=findall(r"[0-9\.]+",data[0])[0]
    if version!="1.2":
        print("unknown file or version (only support v1.2)")
        return False
    
    # decide position of each tag
    tagLattice=tagSupercell=tagOrbitals=tagBonds=tagTemperature=tagSweeps=tagModel=tagAlgorithm=tagNcores=0
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
        if keyword[0]=='Sweeps':
            tagSweeps=iline
        if keyword[0]=='Model':
            tagModel=iline
        if keyword[0]=='Algorithm':
            tagAlgorithm=iline
        if keyword[0]=='Mesurement':
            tagMesurement=iline
        if keyword[0]=='Ncores':
            tagNcores=iline

    if tagLattice*tagSupercell*tagOrbitals*tagBonds*tagTemperature*tagSweeps*tagModel*tagAlgorithm*tagNcores==0:
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
                         [float(ele[1]),float(ele[2]),float(ele[3])],
                         [int(ele[4]),int(ele[5]),(float(ele[6]),float(ele[7]),float(ele[8]))]
                        ])
                         # source    target      [overlat.                                 ] Jz            Jx            Jy                           
        bondList.append([int(ele[4]),int(ele[5]),[float(ele[6]),float(ele[7]),float(ele[8])],float(ele[1]),float(ele[2]),float(ele[3])])

    # load mesurements
    GcPack=findall(r'[0-9\-]+',data[tagMesurement+1])
    s, t, v1, v2, v3 = [int(x) for x in GcPack]
    GcOrb=[[s,t],[v1,v2,v3]]
    # load other parameters
    Tpack=findall(r"[0-9\.]+",data[tagTemperature+1])
    T0, T1, nT = float(Tpack[0]), float(Tpack[1]), int(Tpack[2])
    nTermSweep=findall(r"[0-9\.\-]+",data[tagSweeps+1])
    nthermal, nsweep, ninterval=[int(x) for x in nTermSweep]
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
    gui.corrGui.setValue(GcPack)
    gui.MCparamGui.setValue(nTermSweep)
    gui.modelStr.set(modelType)
    gui.algoStr.set(algorithm)
    gui.coreGui.setValue(_ncores)
    gui.updateStructureViewer()
    return True
    