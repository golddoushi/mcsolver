from tkinter import Label, LabelFrame, Frame, Spinbox, Button, END, VERTICAL, N, S, W, E, StringVar
from multiprocessing import cpu_count
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg#,NavigationToolbar2Tk
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
import toolbox as toolbox
import WannierKit as wan

global gui  # root gui
global latticeGui, supercellGui, latticeData # read lattice matrix
global nOrbnBondGui, nOrb, nBonds  # read number of orbitals and bonds

global OrbListBox, IDandTypeNote, PosNote
global BondBox, IDandTypeOfBondNote, BondDetailNote, AnisotropyOfBondNote

global TlistGui, MCparamGui, modelGui, algorithmGui, coreGui
global resultViewerBase, resultViewer, structureFrame, structureViewer, structureAxis
global submitBtn

#gui=tk.Tk(className='mc solver v0.0.1')

###################
# latice settings #
###################

def loadLatticePannel():
    global gui, latticeGui, supercellGui
    LatticeFrame=LabelFrame(gui,text='Lattice')
    LatticeFrame.grid(row=0,column=0)

    a0_base=Frame(LatticeFrame)
    noteFrame0=toolbox.NoteFrm(a0_base, init_notes=['a1:','',''],init_data=[1,0,0],row=True)
    a0_base.grid(row=0,column=0)

    a1_base=Frame(LatticeFrame)
    noteFrame1=toolbox.NoteFrm(a1_base, init_notes=['a2:','',''],init_data=[0,1,0],row=True)
    a1_base.grid(row=1,column=0)

    a2_base=Frame(LatticeFrame)
    noteFrame2=toolbox.NoteFrm(a2_base, init_notes=['a3:','',''],init_data=[0,0,1],row=True)
    a2_base.grid(row=2,column=0)

    latticeGui=[noteFrame0,noteFrame1,noteFrame2]

    supercell_base=Frame(LatticeFrame)
    supercellGui=toolbox.NoteFrm(supercell_base,init_notes=['SC:','x','x'],init_data=[16,16,1],row=True,entryWidth=3)
    supercell_base.grid(row=0,column=1,sticky='SE')
    

def updateLatticeData():
    global latticeGui, latticeData
    latticeData=[]
    for reporter in latticeGui:
        latticeData.append(reporter.report())
    latticeData=np.array(latticeData)

####################
# orbital settings #
####################

def correspondToOrbList(*arg):
    global OrbListBox, IDandTypeNote, PosNote
    data=OrbListBox.report()
    if len(data)>0:
        IDandTypeNote.entry_list[0].config(state='normal')
        IDandTypeNote.setValue([data[0],data[1],data[2]]) # id, type, spin
        IDandTypeNote.entry_list[0].config(state='disabled')
        PosNote.setValue(data[3])                         # fractional coordinates
        updateStructureViewer(lightOrb=data[0])

def addOrb():
    global OrbListBox, IDandTypeNote, PosNote
    newData=list(OrbListBox.infoData)
    idAndType=IDandTypeNote.report()
    pos=PosNote.report()
    #idAndType.append(pos)
    newData.append([len(newData),int(idAndType[1]),idAndType[2],tuple(pos)])
    OrbListBox.updateInfo(newData)
    updateStructureViewer()
    
def deletOrb():
    global OrbListBox, IDandTypeNote, PosNote
    newData=list(OrbListBox.infoData)
    if len(newData)==0:
        return
    idAndType=IDandTypeNote.report()
    idxs=int(idAndType[0])
    newData.pop(idxs)
    OrbListBox.updateInfo(newData)
    updateStructureViewer()

def resetOrb():
    global OrbListBox, IDandTypeNote, PosNote
    newData=list(OrbListBox.infoData)
    if len(newData)==0:
        return
    idAndType=IDandTypeNote.report()
    pos=PosNote.report()
    idxs=int(idAndType[0])
    #print(idxs)
    newData.pop(idxs)
    newData.insert(idxs,[int(idAndType[0]),int(idAndType[1]),idAndType[2],tuple(pos)])
    OrbListBox.updateInfo(newData)
    updateStructureViewer()

def orbitalDataFormat(info):
    return 'ID: %d Type: %d Spin: %.1f FracX: %.4f %.4f %.4f'%(info[0],info[1],info[2],info[3][0],info[3][1],info[3][2])

def loadOrbitals():
    global gui, OrbListBox, IDandTypeNote, PosNote
    OrbFrame=LabelFrame(gui,text='Orbital list')
    OrbFrame.grid(row=1,column=0,columnspan=1)

    list_base=Frame(OrbFrame)
    list_base.grid(row=0,column=0,columnspan=2)
    OrbListBox=toolbox.InfoList(list_base, correspondToOrbList, orbitalDataFormat, initialInfo=[[0,0,1,(0.,0.,0.)]],width=45,height=5)

    addOrbFrameBase=Frame(OrbFrame)
    addOrbFrameBase.grid(row=1,column=0)

    id_base=Frame(addOrbFrameBase)
    id_base.grid(row=0,column=0)
    IDandTypeNote=toolbox.NoteFrm(id_base, init_notes=['ID:','Type:','Init spin:'],init_data=[0,0,1],row=True,entryWidth=5)
    IDandTypeNote.entry_list[0].config(state='disabled')
    pos_base=Frame(addOrbFrameBase)
    pos_base.grid(row=1,column=0,sticky='W')
    PosNote=toolbox.NoteFrm(pos_base, init_notes=['pos','',''],init_data=[0.,0.,0.],row=True,entryWidth=6)

    addBtn=Button(addOrbFrameBase,text='add',command=addOrb)
    addBtn.grid(row=0,column=1,rowspan=2)
    resetBtn=Button(addOrbFrameBase,text='reset',command=resetOrb)
    resetBtn.grid(row=0,column=2,rowspan=2)
    delBtn=Button(addOrbFrameBase,text='delet',command=deletOrb)
    delBtn.grid(row=0,column=3,rowspan=2)


##################
# Bonds settings #
##################

def correspondToBondList(*arg):
    global BondBox, IDandTypeOfBondNote, BondDetailNote
    data=BondBox.report()
    if len(data)>0:
        IDandTypeOfBondNote.entry_list[0].config(state='normal')
        IDandTypeOfBondNote.setValue([data[0],data[1][0],data[1][1],data[1][2]]) # id, [Jz, Jx, Jy]
        IDandTypeOfBondNote.entry_list[0].config(state='disabled')
        BondDetailNote.setValue([data[2][0],data[2][1],data[2][2][0],data[2][2][1],data[2][2][2]])        # fractional coordinates
        updateStructureViewer(lightID=data[0])

def addBond():
    global BondBox, IDandTypeOfBondNote, BondDetailNote
    newData=list(BondBox.infoData)
    idAndType=IDandTypeOfBondNote.report()
    bondDetail=BondDetailNote.report()
    newData.append([len(newData),
                    [idAndType[1],idAndType[2],idAndType[3]],
                    [int(bondDetail[0]),int(bondDetail[1]),[int(bondDetail[2]),int(bondDetail[3]),int(bondDetail[4])]]])
    BondBox.updateInfo(newData)
    updateStructureViewer()
    
def deletBond():
    global BondBox, IDandTypeOfBondNote, BondDetailNote
    newData=list(BondBox.infoData)
    if len(newData)==0:
        return
    idAndType=IDandTypeOfBondNote.report()
    idxs=int(idAndType[0])
    newData.pop(idxs)
    BondBox.updateInfo(newData)
    updateStructureViewer()

def resetBond():
    global BondBox, IDandTypeOfBondNote, BondDetailNote
    newData=list(BondBox.infoData)
    if len(newData)==0:
        return
    idAndType=IDandTypeOfBondNote.report()
    bondDetail=BondDetailNote.report()
    idxs=int(idAndType[0])
    newData.pop(idxs)
    newData.insert(idxs,[len(newData),
                    [idAndType[1],idAndType[2],idAndType[3]],
                    [int(bondDetail[0]),int(bondDetail[1]),[int(bondDetail[2]),int(bondDetail[3]),int(bondDetail[4])]]])
    BondBox.updateInfo(newData)
    updateStructureViewer()

def bondDataFormat(info):
    return 'ID: %d J: %.3f source: %d target: %d overLat: %d %d %d'%(info[0],info[1][0],info[2][0],info[2][1],info[2][2][0],info[2][2][1],info[2][2][2])

def loadBonds():
    global gui, BondBox, IDandTypeOfBondNote, BondDetailNote, AnisotropyOfBondNote
    BondFrame=LabelFrame(gui,text='Bond list')
    BondFrame.grid(row=2,column=0,columnspan=1)

    list_base=Frame(BondFrame)
    list_base.grid(row=0,column=0,columnspan=2)
    BondBox=toolbox.InfoList(list_base, correspondToBondList, bondDataFormat, 
                             initialInfo=[[0,[-1,-1,-1],[0,0,(1,0,0)]],[1,[-1,-1,-1],[0,0,(0,1,0)]]],
                             width=45,height=5)

    addBondFrameBase=Frame(BondFrame)
    addBondFrameBase.grid(row=1,column=0)

    id_base=Frame(addBondFrameBase)
    id_base.grid(row=0,column=0)
    IDandTypeOfBondNote=toolbox.NoteFrm(id_base, init_notes=['ID:','Jz','Jx','Jy'],init_data=[1,-1,-1,-1],row=True,entryWidth=5)
    IDandTypeOfBondNote.entry_list[0].config(state='disabled')

    detail_base=Frame(addBondFrameBase)
    detail_base.grid(row=1,column=0,sticky=(W,E))
    BondDetailNote=toolbox.NoteFrm(detail_base, init_notes=['s','t','over lat.','',''],init_data=[0,0,1,0,0],row=True,entryWidth=2)

    anis_base=Frame(addBondFrameBase)
    anis_base.grid(row=2,column=0)
    AnisotropyOfBondNote=toolbox.NoteFrm(anis_base, init_notes=['Orb Ani:','Dz','Dx','Dy'],init_data=[0,0,0,0],row=True,entryWidth=4)

    addBtn=Button(addBondFrameBase,text='add',command=addBond)
    addBtn.grid(row=0,column=1,rowspan=2,sticky='E')
    resetBtn=Button(addBondFrameBase,text='reset',command=resetBond)
    resetBtn.grid(row=0,column=2,rowspan=2,sticky='E')
    delBtn=Button(addBondFrameBase,text='delet',command=deletBond)
    delBtn.grid(row=0,column=3,rowspan=2,sticky='E')

###############
# MC settings #
###############

def loadMCSettings():
    global gui, TListGui, MCparamGui, modelGui, algorithmGui, coreGui
    SettingFrame=LabelFrame(gui,text='Other settings')
    SettingFrame.grid(row=3,column=0,sticky=(W,E))

    temp_base=Frame(SettingFrame)
    temp_base.grid(row=0,column=0)
    TListGui=toolbox.NoteFrm(temp_base, init_notes=['T start:','T end','total points:'], init_data=[2.0,2.4,20],row=True,entryWidth=6)

    MCparam_base=Frame(SettingFrame)
    MCparam_base.grid(row=1,column=0,sticky='W')
    MCparamGui=toolbox.NoteFrm(MCparam_base, init_notes=['nthermal:','nsweep:'], init_data=[20000,40000],row=True)

    model_base=Frame(SettingFrame)
    model_base.grid(row=2,column=0,sticky='W')
    label1=Label(model_base,text='Model:')
    label1.grid(row=0,column=0)
    modelGui=Spinbox(model_base,from_=1, to=3, values=['Ising','XY','Heisenberg'],width=12)
    modelGui.grid(row=0,column=1)
    
    label2=Label(model_base,text='Algorithm:')
    label2.grid(row=0,column=2)
    algorithmGui=Spinbox(model_base,from_=1, to=3, values=['Wolff','Metroplis','Sweden-Wang'],width=12)
    algorithmGui.grid(row=0,column=3)

    core_base=Frame(SettingFrame)
    core_base.grid(row=3,column=0,sticky='W')
    coreGui=toolbox.NoteFrm(core_base, init_notes=['core:'], init_data=[np.max([1,int(cpu_count()/2)])])

########################
# Structure visualizer #
########################

def loadStructureViewer():
    global gui, latticeGui, OrbListBox, BondBox, supercellGui, structureFrame, structureViewer, structureAxis
    structureFrame=LabelFrame(gui, text='Structure viewer')
    structureFrame.grid(row=0,column=1,rowspan=2)

    tb=wan.TBmodel()
    a1=latticeGui[0].report()
    a2=latticeGui[1].report()
    a3=latticeGui[2].report()
    tb.lattice=np.array([a1,a2,a3])
    tb.orbital_coor=[[np.array(ele[3]),50,'red'] for ele in OrbListBox.infoData]
    tb.norbital=len(tb.orbital_coor)
    tb.hopping=[[bond_data[2][0],bond_data[2][1],np.array(bond_data[2][2]),bond_data[1][0],'green', 2] for bond_data in BondBox.infoData]
    tb.nhoppings=len(tb.hopping)
    Lx, Ly, Lz=[4 if x>1 else 1 for x in supercellGui.report() ]
    tb.make_supercell([Lx,0,0],[0,Ly,0],[0,0,Lz])
    f, structureAxis=tb.viewStructure()
    
    structureViewer=FigureCanvasTkAgg(f,structureFrame)
    structureViewer.draw()
    structureViewer.get_tk_widget().grid(row=0,column=0)
    structureAxis.figure.canvas=structureViewer
    structureAxis.mouse_init()

def updateStructureViewer(lightID=-1,lightOrb=-1):
    global gui, latticeGui, OrbListBox, BondBox, supercellGui, structureFrame, structureViewer, structureAxis
    elev, azim = structureAxis.elev, structureAxis.azim

    tb=wan.TBmodel()
    a1=latticeGui[0].report()
    a2=latticeGui[1].report()
    a3=latticeGui[2].report()
    tb.lattice=np.array([a1,a2,a3])
    tb.orbital_coor=[[np.array(ele[3]),100 if ele[0]==lightOrb else 50,'yellow' if ele[0]==lightOrb else 'red' ] for ele in OrbListBox.infoData]
    tb.norbital=len(tb.orbital_coor)
    tb.hopping=[[bond_data[2][0],bond_data[2][1],np.array(bond_data[2][2]),bond_data[1][0],'yellow' if bond_data[0]==lightID else 'green', 6 if bond_data[0]==lightID else 2] for bond_data in BondBox.infoData]
    tb.nhoppings=len(tb.hopping)
    Lx, Ly, Lz=[4 if x>1 else 1 for x in supercellGui.report() ]
    tb.make_supercell([Lx,0,0],[0,Ly,0],[0,0,Lz])
    f, structureAxis=tb.viewStructure()
    structureAxis.view_init(elev=elev,azim=azim)
    
    structureViewer.get_tk_widget().destroy()
    structureViewer=FigureCanvasTkAgg(f,structureFrame)
    structureViewer.draw()
    structureViewer.get_tk_widget().pack()
    structureAxis.figure.canvas=structureViewer
    structureAxis.mouse_init()

#####################
# Resutl visualizer #
#####################

def loadResultViewer():
    global gui, resultViewerBase, resultViewer
    resultViewerBase=LabelFrame(gui, text='Result viewer')
    resultViewerBase.grid(row=2,column=1,rowspan=2)

    f=Figure(figsize=(4,3))
    f.add_subplot(111).plot([0,2],[0,0],color='black')
    resultViewer=FigureCanvasTkAgg(f,resultViewerBase)
    resultViewer.draw()
    resultViewer.get_tk_widget().pack()

def updateResultViewer(TList=[],magList=[]):
    global gui, resultViewerBase, resultViewer
    #print('updating:',TList,magList)
    f=Figure(figsize=(4,3))
    ax=f.add_subplot(111)
    ax.scatter(TList,magList,color='red')
    ax.set_title('<Spin> per cell')
    resultViewer.get_tk_widget().destroy()
    resultViewer=FigureCanvasTkAgg(f,resultViewerBase)
    resultViewer.draw()
    resultViewer.get_tk_widget().pack()

#############
# start btn #
#############

def loadStartBtn(submitFunc):
    global gui, submitBtn
    submit_base=Frame(gui)
    submit_base.grid(row=4,column=0,columnspan=2)

    note=Label(submit_base, text='Thanks for your attention and wish you find sth. helpful. \n\
                        Please cite Ref: Magnetic switches via electric field in BN nanoribbons. Applied Surface Science 480(2019)\n\
                        Thank you very much!')
    note.grid(row=0,column=1)
    submitBtn=Button(submit_base,text='Submit',command=submitFunc)
    submitBtn.grid(row=0,column=0)

def loadEverything(root,submitFunc):
    global gui
    gui=root
    loadLatticePannel()
    updateLatticeData()

    loadOrbitals()
    loadBonds()
    loadMCSettings()
    loadStructureViewer()
    loadResultViewer()
    loadStartBtn(submitFunc)
