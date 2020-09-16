'''
Created on 2019 7 22

@author: Andrew
'''
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np 
#import re
import auxiliary as aux
#import paramReader as param
#import posReader as pR
#import dftMain as dft
#import subprocess
import WannierKit
import scipy.optimize as opt
#import matplotlib
#matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt
import fileio as io

global tb       # tight-binding model (for spin-wave calc.)
global S        # local spin
global Jz_eff, Jxy, A_eff  # the coupling constants extracted from DFT calc.
global hBZ      # half the Brillouin zone
global dhBZ     # double the BZ
global norb     # number of mag. ions in cell
global eig0_set # the eigenvalues on the hBZ
global magList_HF  # store the MT curv calc. by Hatree-Fock
global magList_MF  #           MT            by mean-field
global eig0_HF, jka, bka, twojka, JB0ab, B0, JBkaa, fourA, twoA # Hatree-Fock needed
global Jekn, onsite, Jkaa, JBkaa_db
global x0
global Tc, Tc_MF

def mainLoop(rpath):
    global tb,eig0_set,Tc,Tc_MF,magList_HF,magList_MF,S
    # create mc main directory
    print('start renormalized spin-wave theoretical calculations')
    io.loadParam(updateGUI=False,rpath=rpath)
    # initialize tight-binding model
    tb=WannierKit.TBmodel()
    __tbInit()
    tb.plotbands()
    eig0_min=np.min(eig0_set)
    print('spin-wave gap at 0K is: %.6f'%eig0_min)
    if(eig0_min>0):
        Sigma=0
        for eig in eig0_set:
            Sigma+=1./eig
        Tc_=(S+1)*len(eig0_set)/3./Sigma
        print('mean field estimation for Tc: %.3f'%Tc_)
        exit()

        print('start scf procedure for each temperature, to get M-T curv')
        # calc M-T curv
        Tc, magList_HF = __MTCurv(path='./',draw=True,algo='HF_longWave')
        print('Hatree-Fock and long wave approximation give Tc:',Tc,'(Kelvin)')
        print
        Tc_MF,magList_MF=__MTCurv(path='./',draw=True,algo='meanfield')
        print('Mean-field approximation give Tc:',Tc_MF,'(Kelvin)')
        print('Spin-wave theory Tc (Hatree-Fock long wave limit): {:.1f}\n'.format(Tc))
        print('Spin-wave theory Tc (Mean-field approximation): {:.1f}\n'.format(Tc_MF))
        __spinWave_T_vs_Occ('./',algo='HF_longWave')
        __spinWave_T_vs_Occ('./',algo='Mean_Field')
        return Tc
    else:
        print('find zero or negative spin-wave gap!')
        print('Spin-wave theory find no Tc!\n')
        return 0

def __tbInit():
    global tb, hBZ, dhBZ, S, norb
    # parameters for Hatree-Fock calc.
    global eig0_set, fourA, twoA, A

    global eig0_HF, Jekn, onsite, Jkaa, JBkaa_db, x0
    

    # set lattice
    tb.lattice=np.array(io.LMatrix)
    # set kpath for test
    tb.genReciLattice()
    tb.autoGenerateKpath2D(20)
    # get coupling constants
    #Jz, Jxy, A = __genXXYParam()
    #J=Jxy
    #B=Jz-Jxy
    S=io.S[0]
    A=io.DList[0][0]-(io.DList[0][1]+io.DList[0][2])/2
    #print('spin-wave model parameters:')
    #print('J:', J*1000)
    #print('B:', B*1000)
    #print('A', A*1000)
    # set tb hoppings and diagonal terms according to eq. (S9)
    tb.norbital=len(io.S)
    #print(io.S)
    norb=tb.norbital
    onsite=np.ones(tb.norbital)*2*A
    #print(onsite)
    hopping=[]

    Jrab, Brab=[], []
    for bond in io.bondList:
        sourceID, targetID=bond[0], bond[1]
        overLat=bond[2]
        Jz, Jx, Jy= bond[3], bond[4], bond[5]
        J=(Jx+Jy)/2.
        B=Jz-J
        # onsite energy
        onsite[sourceID]+=Jz
        onsite[targetID]+=Jz
        # hopping 
        hopping.append([sourceID,targetID,np.array(overLat),J])
        Jrab.append([sourceID,targetID,np.array(overLat),J])
        Brab.append([sourceID,targetID,np.array(overLat),B])

    # set half Brillouin zone
    hBZ=[]
    dk=1./io.LPack[0]
    for ikptx in range(io.LPack[0]):
        for ikpty in range(int(io.LPack[0]/2)):
            hBZ.append([ikptx*dk,ikpty*dk])
    dhBZ=list(hBZ) 
    dhBZ.extend(hBZ)

    # get Jekn
    tb.onsite_energy=np.zeros(tb.norbital)
    tb.hopping=hopping
    tb.fixHopping()
    tb.constructHam()
    Jekn=np.array([tb.solveHk(kpt=kpt) for kpt in hBZ])

    # jkaa in hBZ
    Jkaa=np.array([np.diag(tb.constructHk(kpt))[0] for kpt in hBZ]).real

    # produce eig0_HF
    eig0_HF=S*(Jekn-onsite)

    # j(k1-k,aa) in double hBZ
    Jkaa_db=np.array([np.diag(tb.constructHk(kpt))[0] for kpt in dhBZ])

    # B(k1-k,aa) in double hBZ
    tb.hopping=Brab
    tb.fixHopping()
    tb.constructHam()
    Bkaa_db=[np.diag(tb.constructHk(kpt))[0] for kpt in dhBZ]

    # J + B
    JBkaa_db=(Jkaa_db+Bkaa_db).real

    # x0 initial occ.
    x0=np.array([0.]*len(eig0_HF))

    fourA=4*A
    twoA=2*A

    # get eig0_set
    #print(onsite)
    #exit()
    tb.onsite_energy=-S*onsite/11.58875
    tb.hopping=hopping
    for hop in tb.hopping:
        hop[3]*=S/11.58875
    tb.fixHopping()
    tb.constructHam()

    __getEigen0()
    
def __getEigen0():
    global tb, hBZ, eig0_set, eig0_HF
    #eig0_HF=np.array([tb.solveHk(kpt=kpt) for kpt in hBZ])
    eig0_set=np.array(list(eig0_HF)).flatten()

def __spinWaveSCF(T=1.):
    '''mean field approximation'''
    global S, eig0_set
    beta=1./T
    def bzSum(eigList,beta):
        ntot=0.
        for eig in eigList:
            ntot+=1./(np.exp(beta*eig)-1.)
        return ntot
    
    N=len(eig0_set)
    n_avg=0
    
    def scf_meanfield(n_avg):
        eigList=(1-n_avg/2/S)*eig0_set
        ntot=bzSum(eigList,beta)
        return (ntot/N-n_avg)**2
    
    out=opt.minimize_scalar(scf_meanfield,bounds=(0,S))
    return out['fun'], out['x']

def __spinWaveSCF_HF(T=1.):
    '''Hatree-Fock level discussed in H1 expression'''
    global S, eig0_HF, Jkaa, JBkaa_db, Jekn, onsite, norb, A, twoA, fourA
    global x0
    beta=1./T
    N=len(eig0_HF)

    def scf_Hatree_Fock(nk_avg):
        ntot=np.sum(nk_avg)
        # independent contribution
        contr_ind=A+ntot/N*(-Jekn+onsite+twoA)  # 
        # dependent contribution
        contr_dep=1./N*np.array([np.dot((JBkaa_db[N-ik:2*N-ik]-Jkaa),nk_avg) for ik in range(N)])
        # renormalized eig
        eig_ren=(eig0_HF+contr_ind).T+contr_dep

        # occ renorm
        nk_renorm=[np.sum(nk)/norb for nk in 1./(np.exp(beta*eig_ren.T)-1.)]
        diff_nk=np.array(nk_renorm)-nk_avg
        return np.dot(diff_nk,diff_nk)

    out=opt.minimize(scf_Hatree_Fock,x0,bounds=[(0,None) for i in range(N)])
    # update x0
    x0=out['x']
    return out['fun'], np.sum(out['x'])*norb/N

def __spinWaveSCF_HF_longWave(T=1.):
    '''Hatree-Fock level using long-wave approximation'''
    global S, eig0_HF, Jkaa, JBkaa_db, Jekn, onsite, norb, A, twoA, fourA
    beta=1./T
    N=len(eig0_HF)

    def scf_HF_longwave(nk_avg):
        # independent contribution
        contr_ind=A+nk_avg*(-Jekn+onsite+twoA)  # modified ignore the A term
        # dependent contribution
        contr_dep=nk_avg*(JBkaa_db[0:N]-Jkaa[0])
        # renormalized eig
        eig_ren=(eig0_HF+contr_ind).T+contr_dep

        # occ renorm
        nk_renorm=np.sum([np.sum(nk) for nk in 1./(np.exp(beta*eig_ren.T)-1.)])/norb/N
        diff_nk=nk_renorm-nk_avg
        return diff_nk*diff_nk

    out=opt.minimize_scalar(scf_HF_longwave,bounds=(0,S/norb),method='Bounded')
    return out['fun'], out['x']*norb

def __spinWave_T_vs_Occ(path,algo='HF_longWave'):
    '''Draw the T vs Occ. phase'''
    global S, Tc, Tc_MF, eig0_HF, eig0_set, Jkaa, JBkaa_db, Jekn, onsite, norb, A, twoA, fourA
    N=len(eig0_HF)

    def scf_HF_longwave(nk_avg,beta):
        # independent contribution
        contr_ind=A+nk_avg*(-Jekn+onsite+twoA)  # modified ignore the A term
        # dependent contribution
        contr_dep=nk_avg*(JBkaa_db[0:N]-Jkaa[0])
        # renormalized eig
        eig_ren=(eig0_HF+contr_ind).T+contr_dep

        # occ renorm
        nk_renorm=np.sum([np.sum(nk) for nk in 1./(np.exp(beta*eig_ren.T)-1.)])/norb/N
        diff_nk=nk_renorm-nk_avg
        return np.min([-np.log(np.max([diff_nk*diff_nk,1E-10])),10])

    def scf_meanfield(n_avg,beta):
        eigList=(1-n_avg/2/S)*eig0_set
        ntot=np.sum(1./(np.exp(beta*eigList)-1.))
        diff=ntot/len(eigList)-n_avg
        return np.min([-np.log(np.max([diff*diff,1E-10])),10])

    T=np.linspace(1,Tc*1.2,101)
    T_mf=np.linspace(1,Tc_MF*1.2,101)
    Occ=np.linspace(0,S/norb,101)
    Occ_tot=np.linspace(0,S,101)
    phase=[]
    if algo=='HF_longWave':
        for t in T:
            phase.append([scf_HF_longwave(nk_avg,beta=1./t) for nk_avg in Occ])
        xlabel=[int(num) for num in np.linspace(1,Tc*1.2,6)]
    else:
        for t in T_mf:
            phase.append([scf_meanfield(n_avg,beta=1./t) for n_avg in Occ_tot])
        xlabel=[int(num) for num in np.linspace(1,Tc_MF*1.2,6)]
    plt.figure(figsize=(6,6))
    plt.imshow(np.array(phase).T,origin='low',extent=(0,10,0,10))
    #plt.colorbar()
    ylabel=['%.2f'%num for num in np.linspace(0,S,6)]
    plt.xticks(range(0,12,2),xlabel)
    plt.yticks(range(0,12,2),ylabel)
    #plt.xlabel(r'Temperature (K)')
    #plt.ylabel(r'Spin-wave occ.')
    plt.tick_params(labelsize=20)
    plt.savefig(path+algo+'OT.png',dpi=300)
    plt.close()

def __MTCurv(path='./',draw=False,algo='HF_longWave'):
    global S, eig0_set, norb
    global eig0_HF
    global twoA,fourA,A
    global JBkaa_db, Jekn, Jkaa, onsite 

    magList=[]
    TList=[]
    T=1.
    mag0=-5
    while(True):
        if(algo=='HF_longWave'):
            fun, n_avg = __spinWaveSCF_HF_longWave(T)
        elif(algo=='HF'):
            fun, n_avg = __spinWaveSCF_HF(T)
        else:
            fun, n_avg = __spinWaveSCF(T)
        mag=2*(S-n_avg)
        print(T, fun, mag, mag-mag0)
        mag0=mag
        if(fun>1e-8):
            break
        if(mag<0):
            break
        magList.append(mag)
        TList.append(T)
        T+=1
    if(draw):
        plt.figure()
        plt.scatter(TList,magList)
        plt.xlim(0,T)
        plt.savefig(path+algo+'MT.png',dpi=300)
        plt.close()
    return T, magList

mainLoop('./samples/Fe110')