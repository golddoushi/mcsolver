'''
Created on 2018 12 8 

@author: Andrew
'''
import numpy as np
import re
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.axes3d import Axes3D

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