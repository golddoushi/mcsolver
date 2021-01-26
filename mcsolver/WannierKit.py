'''
Created on 2018 12 8 

@author: Andrew
'''
import numpy as np
import re
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
try:
    from . import auxiliary as aux
except:
    import auxiliary as aux

class TBmodel(object):
    Ham=[]
    Ham_k=[]
    lattice=[]
    reci_lattice=[]
    norbital=0
    orbital_coor=[]
    onsite_energy=[]
    
    nhoppings=0
    hopping=[]
    
    kpath=[]
    kmesh=[]
    
    wcc_path=[]
    wcc_polar_direction=-1
    nelectrons=-1
    
    def __init__(self):
        pass
    
    def make_supercell(self,sc_dir0=[1.,1.,0.],sc_dir1=[1.,-1.,0.],sc_dir2=[0.,0.,1.],toHome=True):
        '''
        construct a supercell model.
        often used with slab fun. to get slab model with chosen border.
        '''
        Tmatrix=np.array([sc_dir0,sc_dir1,sc_dir2])  # sc_lat = T dot lat 
        invTmatrix=np.linalg.inv(Tmatrix)            #    lat = invT dot sc_lat
        # construct lattice matrix for supercell
        sc_lattice=np.dot(Tmatrix,self.lattice)
    
        # construct orbitals for supercell
        R0max=int(np.max(abs(Tmatrix[:,0])))
        R1max=int(np.max(abs(Tmatrix[:,1])))
        R2max=int(np.max(abs(Tmatrix[:,2])))
        
        ### find which sub-cell is interior, just check the origins
        interior_cell_list=[]
        interior_cell_coor=[]
        for a0 in range(-R0max,R0max+1):
            for a1 in range(-R1max,R1max+1):
                for a2 in range(-R2max,R2max+1):
                    ori_lat=np.array([a0,a1,a2])
                    ori_sc_lat=np.dot(ori_lat,invTmatrix)
                    interior=True
                    for ori_coor in ori_sc_lat:
                        if(ori_coor<0 or ori_coor>=1):
                            interior=False
                    if(interior):
                        interior_cell_list.append(ori_lat)
                        interior_cell_coor.append(ori_sc_lat)
        
        ### now add every orbs and their attributes
        sc_orbital_coor=[]
        onsite_energy=[]
        norbital=0
        
        hopping=[]
        nhoppings=0
        for iinter_cell, interior_cell in enumerate(interior_cell_list):
            for iorb0, orb0 in enumerate(self.orbital_coor):
                # add coordinates
                orb_sc=np.dot(orb0[0]+interior_cell,invTmatrix)
                sc_orbital_coor.append([orb_sc,orb0[1],orb0[2]])
                norbital+=1
                # add hopping
                for hopping_ele in self.hopping:
                    if(hopping_ele[0]==iorb0):
                        iorb0_sc=iinter_cell*self.norbital+iorb0
                        #hopping_orb0=[iorb0_sc]
                        iorb1, aug_vec, amplify, color, linewidth= hopping_ele[1], hopping_ele[2], hopping_ele[3], hopping_ele[4], hopping_ele[5]
                        #print('iorb0, iorb1:', iorb0, iorb1, 'aug_vec:',aug_vec)
                        orb1_cell_ori=interior_cell+aug_vec # the origin coordinates (a1,a2,a3 lattice)
                        orb1_cell_ori_sc=np.dot(orb1_cell_ori,invTmatrix) # in supercell lattice
                        orb1_cell_ori_sc_reduced=orb1_cell_ori_sc%1       # reduced by 1
                        aug_vec_sc=orb1_cell_ori_sc-orb1_cell_ori_sc_reduced # aug vector in supercell frame
                        # check which level is reached
                        for iori_coor, ori_coor in enumerate(interior_cell_coor):
                            if(ori_coor==orb1_cell_ori_sc_reduced).all():
                                iorb1_sc=iori_coor*self.norbital+iorb1
                                break
                        #print('SC: iorb0, iorb1:', iorb0_sc, iorb1_sc, 'aug_vec:',aug_vec_sc)
                        hopping_orb0=[iorb0_sc, iorb1_sc, aug_vec_sc, amplify, color, linewidth]
                        hopping.append(hopping_orb0)
                        nhoppings+=1
               
        self.lattice=sc_lattice
        self.orbital_coor=sc_orbital_coor
        self.norbital=norbital
        self.onsite_energy=onsite_energy
        self.nhoppings=nhoppings
        self.hopping=hopping
        if toHome:
            self.toHome()

    def toHome(self):
        '''move all coordinates to the home cell'''
        for iorb, orb in enumerate(self.orbital_coor):
            orb_reduced=orb[0]%1
            orb_shift=orb_reduced-orb[0]
            if (orb_shift).any():
                #print(orb,orb_shift)
                # update orb
                orb[0]%=1
                # update the hoppings
                for hopping in self.hopping:
                    # hopping start with orb
                    if(hopping[0]==iorb):
                        hopping[2]+=orb_shift
                    if(hopping[1]==iorb):
                        hopping[2]-=orb_shift

    def viewStructure(self):
        '''visualize the orbital and bonding's spatial configuration'''
        f=Figure(figsize=(3,3))
        ax=f.add_subplot(111,projection='3d')
        
        def dragLineWithTwoPoints(pt0, pt1, c='black', linewidth=2):
            ax.plot([pt0[0],pt1[0]],[pt0[1],pt1[1]],[pt0[2],pt1[2]],c=c, linewidth=linewidth)
         
        # draw the boder of unit cell   
        pt0=np.zeros(3)
        a1,a2,a3=self.lattice[0],self.lattice[1],self.lattice[2]
        ax.text(a1[0],a1[1],a1[2],'$a_1$')
        ax.text(a2[0],a2[1],a2[2],'$a_2$')
        ax.text(a3[0],a3[1],a3[2],'$a_3$')

        dragLineWithTwoPoints(pt0,a1)
        #dragLineWithTwoPoints(a1,a1+a2)
        #dragLineWithTwoPoints(a1+a2,a2)
        dragLineWithTwoPoints(pt0,a2)
        
        #dragLineWithTwoPoints(a3,a1+a3)
        #dragLineWithTwoPoints(a1+a3,a1+a2+a3)
        #dragLineWithTwoPoints(a1+a2+a3,a2+a3)
        #dragLineWithTwoPoints(a3,a2+a3)
        
        dragLineWithTwoPoints(pt0,a3)
        #dragLineWithTwoPoints(a1,a1+a3)
        #dragLineWithTwoPoints(a2,a2+a3)
        #dragLineWithTwoPoints(a1+a2,a1+a2+a3)
        
        for iorb, orb in enumerate(self.orbital_coor):
            if iorb>=50:
                print('num. of orb >=50, some orbs will not be illustrated')
                break
            orb_xyz=np.dot(orb[0],self.lattice)
            ax.scatter(orb_xyz[0],orb_xyz[1],orb_xyz[2],c=orb[2],s=orb[1])
            #ax.text(orb_xyz[0],orb_xyz[1],orb_xyz[2],str(iorb))
        
        for ihopping, hopping in enumerate(self.hopping):
            if ihopping>=200:
                print('num. of hopping >=200, some hoppings will not be illustrated')
                break
            iorb0, iorb1, Rextra, amp, color, linewidth=hopping[0], hopping[1], hopping[2], hopping[3], hopping[4], hopping[5]
            if np.dot(Rextra,Rextra)>0:
                continue
            orb0_xyz=np.dot(self.orbital_coor[iorb0][0],self.lattice)
            orb1_xyz=np.dot((self.orbital_coor[iorb1][0]+Rextra),self.lattice)
            dragLineWithTwoPoints(orb0_xyz,orb1_xyz,c=color,linewidth=linewidth)
            
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.grid(False)
        ax.axis('off')
        return f, ax

    def genReciLattice(self):
        # generate the 3rd lattice vector if input 2D model
        def antiSymmMatrix(a):
            return np.array([
                [ 0   ,-a[2], a[1]],
                [ a[2], 0   ,-a[0]],
                [-a[1], a[0], 0   ]
                ])
        
        if(len(self.lattice)==2):
            a1=np.concatenate((self.lattice[0],[0.]))
            a2=np.concatenate((self.lattice[1],[0.]))
            a3=np.array([0.,0.,1.])
            self.lattice=np.array([a1,a2,a3])
        a1=self.lattice[0]
        a2=self.lattice[1]
        a3=self.lattice[2]
        
        V=np.abs(np.dot(np.dot(antiSymmMatrix(a1),a2),a3))
        
        b1=np.dot(antiSymmMatrix(a2),a3)
        b2=np.dot(antiSymmMatrix(a3),a1)
        b3=np.dot(antiSymmMatrix(a1),a2)
        
        b1*=2*np.pi/V
        b2*=2*np.pi/V
        b3*=2*np.pi/V
        self.reci_lattice=np.array([b1,b2,b3])

    def fixHopping(self):
        '''
        throw redundant hoppings and add hoppings accordingto time-reversal-symmetry (TRS)
        '''
        perfect_hoppings=[]
        def redundantless(hop_check):
            for hop_ref in perfect_hoppings:
                if(hop_check[0]==hop_ref[0] and hop_check[1]==hop_ref[1] and\
                   (hop_check[2]==hop_ref[2]).all()):
                    if(abs(hop_check[3]-hop_ref[3])>0.01):
                        print ('WARNING: hoppings brreaks time-reversal symmetry!')
                    return False
            return True
        for hopping in self.hopping:
            if redundantless(hopping):
                perfect_hoppings.append(hopping)
                # TRS
                hasTRS=False
                if(hopping[0]!=hopping[1]):
                    hasTRS=True
                elif(np.dot(hopping[2],hopping[2])>0.01):
                    hasTRS=True
                if(hasTRS):
                    tr_hopping=[hopping[1],hopping[0],-hopping[2],hopping[3]]
                    if redundantless(tr_hopping):
                        perfect_hoppings.append(tr_hopping)
        self.hopping=perfect_hoppings

    def constructHam(self):
        '''
        construct real-space Hamiltonian'''
        # firstly, construct the architecture
        self.Ham=[]
        for origin_orbital_index in range(self.norbital):
            hopping_tmp=[]
            for destiny_orbital_index in range(self.norbital):
                hopping_tmp.append([])
            self.Ham.append(hopping_tmp)
        
        # then import all hoppings
        for hopping in self.hopping:
                    # source     # target            # overlat  # t
            self.Ham[hopping[0]][hopping[1]].append([hopping[2],hopping[3]])
            
    def constructHk(self,kpt=[0.,0.,0.],debug=False):
        "Construct Hamiltonian for a certain k-point with reduced coordinates"
        if len(kpt)==2:
            kpt=list(kpt) 
            kpt.append(0)
        kpt=np.array(kpt)
            
        Ham_k=[]
        for ibloch in range(self.norbital):
            Hk_tmp=[]
            for jbloch in range(self.norbital):
                matrix_element=0.
                for hopping in self.Ham[ibloch][jbloch]:
                    matrix_element+=hopping[1]*np.exp(-2j*np.pi*np.dot(kpt,hopping[0]),dtype=complex)
                
                Hk_tmp.append(matrix_element)
            Ham_k.append(Hk_tmp)
        #if debug:
        #    print(Ham_k)
            #print(self.onsite_energy)
        Ham_k=np.array(Ham_k)+self.onsite_energy*np.eye(self.norbital,dtype=complex)
        
        #if debug:
        #    print(Ham_k)
        #    exit()
        return Ham_k
    
    def solveHk(self,kpt=[0.,0.,0.],return_orb=False,debug=False):
        '''solve the eigen-values for a k-point with reduced coordinates,
           and return the eigen-vectors optionally
        '''
        #if debug:
        #    print(kpt)
        #    print(self.Ham)
        #    exit()
        Ham_k=self.constructHk(kpt=kpt,debug=debug)
        #print(Ham_k)
        eig, vec = np.linalg.eigh(Ham_k,'U')
        #if debug: exit()
        
        if return_orb:
            return eig, vec
        else:
            return eig

    def genKPath(self,highSymK,nikpt):
        '''
        generate kpath between high-symmetry kpoints
        '''
        kpath=[]
        npath=len(highSymK)-1
        for ipath in range(npath):
            dk=(highSymK[ipath+1]-highSymK[ipath])/nikpt[ipath]
            for ikpt in range(nikpt[ipath]):
                kpath.append(ikpt*dk+highSymK[ipath])
        kpath.append(highSymK[npath])
        return kpath
    
    def autoGenerateKpath2D(self,nikpt=20):
        '''
        generate k-path automatically,
        only for 2D case, temporally
        Parameter(s):
            nikpt: represents the number of interval KPs between two high-symmetric KP
        '''
        b1=self.reci_lattice[0][:2]/2
        b2=self.reci_lattice[1][:2]/2
        b1b2=b1+b2
        # calc. coordinates of edge
        # in cartesian coor.
        factor=np.array([b1,b2])
        res=np.array([np.dot(b1,b1),np.dot(b2,b2)])
        K_cart=np.dot(np.linalg.inv(factor),res)
        # in reci-fractional coor.
        K0=np.dot(np.array([K_cart[0],K_cart[1],0.]),np.linalg.inv(self.reci_lattice))
        
        # another possible
        factor=np.array([b1,b1b2])
        res=np.array([np.dot(b1,b1),np.dot(b1b2,b1b2)])
        K_cart=np.dot(np.linalg.inv(factor),res)
        K1=np.dot(np.array([K_cart[0],K_cart[1],0.]),np.linalg.inv(self.reci_lattice))
        
        # chose the short one
        K=K0[:2] if np.dot(K0,K0)<np.dot(K1,K1) else K1[:2]
        
        # find the coordinates of bond
        # four possible positions
        M_list=[np.array([0.5,0]),np.array([-0.5,0]),np.array([0.,0.5]),np.array([0.,-0.5])]
        sort_list=aux.quicksort([np.dot(M-K,M-K) for M in M_list])
        # chose the nearest one
        M=M_list[sort_list[0]]
        G=np.array([0.,0.])
        
        # generate kpath
        print('auto. generated high symmetric point in k-space')
        print(G,K,M,G)
        self.kpath=self.genKPath([G,K,M,G],[nikpt,nikpt,nikpt])

    def plotbands(self,nfermi=-1,path='./'):
        "plot electronic band structures"
        "if nfermi>0 then only plot bands above and underneath half-filling for nfermi"
        "e.g., nfermi=1 then only the top valence band and lowest conducting band are plotted"
        self.constructHam()
        plt.figure()
        kpath=np.linspace(0,1,len(self.kpath))
        eig=[]
        for kpt in self.kpath:
            eig_list = self.solveHk(kpt=kpt,debug=True)
            eig.append(eig_list)
            #for eig in eig_list:
            #    plt.scatter(ikpt,eig,color='black')
        eig=np.array(eig)
        if(nfermi==-1):
            for iband in range(self.norbital):
                plt.plot(kpath,eig[:,iband],color='black')
        else:
            lower=int(self.norbital/2-nfermi)
            higher=int(self.norbital/2+nfermi)
            for iband in range(lower,higher):
                plt.plot(kpath,eig[:,iband],color='black')
        plt.xlim(0,1)
        #plt.ylim(0,1.1*np.max(eig.flatten()))
        plt.xticks([])
        plt.savefig(path+'spectra.png',dpi=300)
        plt.close()

    def plot2DStructure(self,vec,Lx=1,Ly=1,kpt=np.array([0,0,0]),S=1.5):
        #print(self.orbital_coor)
        #print(self.lattice)
        #plt.figure()
        fout=open('./spinWave.txt','w')
        fout.write('#X     Y     Z     dX    dY    dZ\n')
        for x in range(Lx):
            for y in range(Ly):
                R=np.array([x,y,0])
                blochPhase=2*np.pi*((kpt.dot(R))%1)
                #print('%d %d %.3f'%(x,y,blochPhase))
                for iorb, orb in enumerate(self.orbital_coor):
                    subOrbPos=(orb[0]+R).dot(self.lattice)

                    subOrbPhase=np.log(vec[iorb]/abs(vec[iorb])).imag if abs(vec[iorb])>1e-5 else 0
                    #phi=(subOrbPhase+blochPhase)#*np.pi*2
                    wycoffPhase=orb[0].dot(kpt)*2*np.pi
                    phi=subOrbPhase+blochPhase+wycoffPhase
                    theta=abs(vec[iorb])*0.5/S*np.pi*0.5
                    nS=(np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta))
                    #print('theta %.3f phi %.3f'%(theta,phi))
                    #print(vec)
                    fout.write('%.6f %.6f %.6f %.6f %.6f %.6f\n'%(*subOrbPos,*nS))
                    #fout.write('%.6f %.6f %.6f %.6f\n'%(*subOrbPos,phi))
                    #plt.scatter(subOrbPos[0],subOrbPos[1],c='black')
                    #plt.annotate('%.3f'%(totalPhase/np.pi/2),(subOrbPos[0],subOrbPos[1]))
                #exit()
        #plt.show()
        fout.close()
