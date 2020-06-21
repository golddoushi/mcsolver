#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Vec
{
    int dimension;
    double x,y,z;
}Vec;
// vec1=abs(vec1)
void normalize(Vec *vec1){
    double len=sqrt(vec1->x*vec1->x+vec1->y*vec1->y+vec1->z*vec1->z);
    if (len<1e-5) return;
    vec1->x/=len;
    vec1->y/=len;
    vec1->z/=len;
}
// vec1=vec2
void equal(Vec *vec1, Vec vec2){
    vec1->x=vec2.x;
    vec1->y=vec2.y;
    vec1->z=vec2.z;
}
// vec1=fabs(vec1)
void vabs(Vec *vec1){
    vec1->x=fabs(vec1->x);
    vec1->y=fabs(vec1->y);
    vec1->z=fabs(vec1->z);
}
// vec1/=c
void cDivides(Vec *vec1, double c){
    vec1->x/=c;
    vec1->y/=c;
    vec1->z/=c;
}
// vec1*=c
void cTimes(Vec *vec1, double c){
    vec1->x*=c;
    vec1->y*=c;
    vec1->z*=c;
}
// vec1+=vec2
void plusEqual(Vec *vec1, Vec vec2){
    vec1->x+=vec2.x;
    vec1->y+=vec2.y;
    vec1->z+=vec2.z;
}
// vec1-=vec2
void minusEqual(Vec *vec1, Vec vec2){
    vec1->x-=vec2.x;
    vec1->y-=vec2.y;
    vec1->z-=vec2.z;
}
// vec1 dot vec2
double dot(Vec vec1, Vec vec2){
    return vec1.x*vec2.x+vec1.y*vec2.y+vec1.z*vec2.z;
}

double diagonalDot(Vec vec1, Vec vec2, Vec vec3){
    return vec1.x*vec2.x*vec3.x+
           vec1.y*vec2.y*vec3.y+
           vec1.z*vec2.z*vec3.z;
}
Vec *generateRandomVec(void){
    Vec *direction=(Vec*)malloc(sizeof(Vec));
    direction->x=rand()/(double) RAND_MAX-0.5;
    direction->y=rand()/(double) RAND_MAX-0.5;
    direction->z=rand()/(double) RAND_MAX-0.5;
    double len2=dot(*direction, *direction);
    if (len2>0.5)
    {
        free(direction);
        return generateRandomVec();
    }else
    {
        double len=sqrt(len2);
        direction->x/=len;
        direction->y/=len;
        direction->z/=len;
        return direction;
    }
}

typedef struct Orb
{
    int id;
    double S; 
    Vec spin;
    Vec transSpin;
    int nlink;
    Vec *linkStrength;
    int inBlock;
    struct Orb **linkedOrb;
    struct Orb **linkedOrb_rnorm;
    Vec onsiteAnisotropy;

    double d_onsiteEnergy;
    double sDotN;
    int isProjected;
    int chosen;
    int nOrbInCluster;
    struct Orb **orb_cluster;
    
    double h;
}Orb;

//establishLattice(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength);
void establishLattice(Orb *lattice, int totOrbs, double initSpin[totOrbs], double initD[totOrbs][3], double h, double flunc, int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking][3],
                      int totOrb_rnorm, int nOrbInCluster, int rOrb[totOrb_rnorm], int rOrbCluster[totOrb_rnorm][nOrbInCluster]){
    //printf("establishing whole lattice with %d orbs and %d linkings for each orb\n",totOrbs,maxNLinking);
    for(int i=0;i<totOrbs;i++){
        //printf("check point-1, entering setting %d",i);
        lattice[i].id=i;
        lattice[i].S=initSpin[i];
        //lattice[i].spin.coor=(double*)malloc(3*sizeof(double));  // allocate spin vector for each orb
        lattice[i].spin.x=initSpin[i];
        lattice[i].spin.y=0;
        lattice[i].spin.z=0;

        Vec *fluncSpin=generateRandomVec();
        cTimes(fluncSpin, flunc);
        plusEqual(&lattice[i].spin,*fluncSpin);
        normalize(&lattice[i].spin);
        cTimes(&lattice[i].spin,fabs(initSpin[i]));
        free(fluncSpin);
        //printf("%.3f %.3f \n",lattice[i].spin.x, fabs(initSpin[i]));
        lattice[i].onsiteAnisotropy.x=initD[i][0];
        lattice[i].onsiteAnisotropy.y=initD[i][1];
        lattice[i].onsiteAnisotropy.z=initD[i][2];

        lattice[i].h=h;
        //printf("orb %d spin: %.3f %.3f %.3f\n",i,lattice[i].spin.coor[0],lattice[i].spin.coor[1],lattice[i].spin.coor[2]);
        lattice[i].transSpin.x=0;  // allocate trans spin vector for each orb
        lattice[i].transSpin.y=0;
        lattice[i].transSpin.z=0;
        lattice[i].nlink=nlink[i];
        //printf("check point 1, orb: %d\n",lattice[i].id);
        lattice[i].linkStrength=(Vec*)malloc(nlink[i]*sizeof(Vec)); // allocate strength for each linking
        //printf("check point 2, total links:%d\n",nlink[i]);
        for(int j=0;j<nlink[i];j++){
            //lattice[i].linkStrength[j].coor=(double*)malloc(3*sizeof(double));
            //printf("check point 3\n");
            lattice[i].linkStrength[j].x=linkStrength[i][j][0];
            lattice[i].linkStrength[j].y=linkStrength[i][j][1];
            lattice[i].linkStrength[j].z=linkStrength[i][j][2];
            //printf("check point 4, link:%d, strength: %.3f %.3f %.3f\n",j,lattice[i].linkStrength[j].coor[0],lattice[i].linkStrength[j].coor[1],lattice[i].linkStrength[j].coor[2]);
        }
        lattice[i].chosen=0;
    }
    //printf("now checking the orb cluster\n");
    for(int i=0;i<totOrb_rnorm;i++){
        int id=rOrb[i];
        lattice[id].chosen=1;
        lattice[id].nOrbInCluster=nOrbInCluster;
        lattice[id].orb_cluster=(Orb**)malloc(nOrbInCluster*sizeof(Orb*));
        //printf("orb%d is chosen to be the center of cluster, involving %d orbs in total:\n",lattice[id].id,lattice[id].nOrbInCluster);
        for(int iorb=0;iorb<nOrbInCluster;iorb++){
            lattice[id].orb_cluster[iorb]=lattice+rOrbCluster[i][iorb];
            //printf("    from input id=%d orb%d\n",rOrbCluster[id][iorb],lattice[id].orb_cluster[iorb]->id);
        }
    }
    //printf("orbitals successfully built\n");
}

void establishLinking(Orb *lattice, int totOrbs, int maxNLinking, int nlink[totOrbs], int linkedOrb[totOrbs][maxNLinking],
                      int totOrb_rnorm, int rOrb[totOrb_rnorm], int linkedOrb_rnorm[totOrb_rnorm][maxNLinking]){
    for(int iorb=0;iorb<totOrbs;iorb++){
        lattice[iorb].linkedOrb=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb][ilink];
        }
    }
    for(int i=0;i<totOrb_rnorm;i++){
        int iorb=rOrb[i];
        lattice[iorb].linkedOrb_rnorm=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb_rnorm[ilink]=lattice+linkedOrb_rnorm[i][ilink];
        }
    }
}

double getCorrEnergy(Orb *source){
    //printf("start calc corr. energy\n");
    double corr=0;
    for(int i=0;i<source->nlink;i++){
        corr+=diagonalDot(source->linkStrength[i],source->spin,source->linkedOrb[i]->spin);
        //printf("J=%.3f S1=%.3f S2=%.3f\n",source->linkStrength[i],source->spin,source->linkedOrb[i]->spin);
        //printf("corr=%.3f\n",corr);
    }
    return corr;
}

double getOnsiteEnergy(Orb *source){  // onsite term
    return diagonalDot(source->onsiteAnisotropy,source->spin,source->spin)+source->spin.z*source->h;
}

Vec getMajoritySpin(Orb*_orb){ // core algorithm for renormalization
    //return _orb->spin; // decimation algorithm
    Vec avgSpin;
    avgSpin.x=0;
    avgSpin.y=0;
    avgSpin.z=0;
    for(int iorb=0;iorb<_orb->nOrbInCluster;iorb++){
        plusEqual(&avgSpin,_orb->orb_cluster[iorb]->spin);
    }
    normalize(&avgSpin);
    cTimes(&avgSpin, _orb->S);
    return avgSpin; // majority algorithm
}

double getCorrEnergy_rnorm(Orb *source){
    double corr=0;
    Vec avgSpin_source=getMajoritySpin(source);
    for(int i=0;i<source->nlink;i++){
        //printf("link to orb%d\n",source->linkedOrb_rnorm[i]->id);
        Vec avgSpin_target=getMajoritySpin(source->linkedOrb_rnorm[i]);
        corr+=diagonalDot(source->linkStrength[i],avgSpin_source,avgSpin_target);
    }
    //printf("Ecorr=%.3f\n",corr);
    return corr;
}

double getOnsiteEnergy_rnorm(Orb *source){
    Vec avgSpin_source=getMajoritySpin(source);
    return source->onsiteAnisotropy.x*avgSpin_source.x*avgSpin_source.x+
           source->onsiteAnisotropy.y*avgSpin_source.y*avgSpin_source.y+
           source->onsiteAnisotropy.z*avgSpin_source.z*avgSpin_source.z+source->h*avgSpin_source.z;
}

double getDeltaCorrEnergy(Orb *source){
    double corr=0;
    for(int i=0;i<source->nlink;i++){
        corr+=diagonalDot(source->linkStrength[i],source->transSpin,source->linkedOrb[i]->spin);
    }
    return corr;
}

double getDeltaOnsiteEnergy(Orb *source){
    double s1x=source->spin.x+source->transSpin.x;
    double s1y=source->spin.y+source->transSpin.y;
    double s1z=source->spin.z+source->transSpin.z;
    //printf("h %.3f, transSpin %.3f\n",source->h,source->transSpin.z);
    return source->onsiteAnisotropy.x*(s1x*s1x-source->spin.x*source->spin.x)+
           source->onsiteAnisotropy.y*(s1y*s1y-source->spin.y*source->spin.y)+
           source->onsiteAnisotropy.z*(s1z*s1z-source->spin.z*source->spin.z)+
           source->h*source->transSpin.z;
}

int expandBlock(int*beginIndex, int*endIndex, Orb *buffer[], int*blockLen, Orb *block[], Vec refDirection){
    //printf("  Buffer: now start and end pt is %d, %d.\n",*beginIndex, *endIndex);
    if(*beginIndex>*endIndex) return 0;

    // FIFO
    Orb *outOrb=buffer[*beginIndex];
    *beginIndex+=1; // pop out the first element
    
    //FILO
    //Orb *outOrb=buffer[*endIndex];
    //*endIndex-=1; // pop out the last element

    if(outOrb->isProjected==0){
        outOrb->sDotN=-2*dot(outOrb->spin,refDirection);
        equal(&outOrb->transSpin, refDirection);
        cTimes(&outOrb->transSpin, outOrb->sDotN);
        outOrb->isProjected=1;
    }
    //printf("center orb: %d\n", outOrb->id);
    //printf("projection along axis: %.3f\n", -s1n/2);
    //printf("variation of spin: %.3f %.3f %.3f\n",outOrb->transSpin.coor[0],outOrb->transSpin.coor[1],outOrb->transSpin.coor[2]);
    //printf("it has %d neighbor orbs\n",outOrb->nlink);
    int i;
    for(i=0;i<outOrb->nlink;i++){
        //double effectiveJ=diagonalDot(outOrb->linkStrength[i], refDirection, refDirection);
        //printf("the %d linking has strength %.3f %.3f %.3f (original) %.3f (effective)\n",i,outOrb->linkStrength[i].coor[0],outOrb->linkStrength[i].coor[1],outOrb->linkStrength[i].coor[2],effectiveJ);
        Orb *linkedOrb=outOrb->linkedOrb[i];
        //printf("consider to add orb %d\n",linkedOrb->id);
        //printf("      considering the %d orb which is linking to %d orb, it is %d in block \n", linkedOrb->id, outOrb->id, linkedOrb->inBlock);
        if(linkedOrb->inBlock==0){
            //printf("projection along axis: %.3f\n", s2n);
            //double corr=-s1n*diagonalDot(refDirection,outOrb->linkStrength[i],linkedOrb->spin); // bond strength
            //printf("      spin of orb %d is %.3f %.3f %.3f and bond strength is %.3f\n",linkedOrb->id,linkedOrb->spin.coor[0],linkedOrb->spin.coor[1],linkedOrb->spin.coor[2],corr);
            if (linkedOrb->isProjected==0)
            {
                linkedOrb->sDotN=-2*dot(linkedOrb->spin,refDirection);
                equal(&linkedOrb->transSpin, refDirection);
                cTimes(&linkedOrb->transSpin, linkedOrb->sDotN);
                linkedOrb->isProjected=1;
            }
            double corr=outOrb->sDotN*linkedOrb->sDotN*diagonalDot(refDirection,refDirection,outOrb->linkStrength[i])/2;
            
            //linkedOrb->d_onsiteEnergy=getDeltaOnsiteEnergy(linkedOrb);
            if(corr<0 && (1-exp(corr))>rand()/(double) RAND_MAX){
                //printf("          -->>fortunately it is added to block with possibility %f\n",(1-exp(2*corr)));
                // update block
                *blockLen+=1;
                block[*blockLen-1]=linkedOrb;
                linkedOrb->inBlock=1;  // register into block
                // update buffer
                *endIndex+=1;
                buffer[*endIndex]=linkedOrb;
            }
        }
    }
    return 1;
}

void blockUpdate(int totOrbs, Orb lattice[], double*p_energy, Vec *p_totSpin){
    //printf("one block update step is initializaing...\n");
    Orb *block[totOrbs];
    Orb *buffer[totOrbs];
    int seedID=rand()%totOrbs;
    block[0]=lattice+seedID;
    buffer[0]=lattice+seedID;
    block[0]->inBlock=1;
    
    int beginIndex=0, endIndex=0, blockLen=1, i, j;
    int *p_beginIndex=&beginIndex, *p_endIndex=&endIndex, *p_blockLen=&blockLen;

    Vec *refDirection=generateRandomVec();
    //refDirection.coor=(double*)malloc(3*sizeof(double));
    //equal(&refDirection, &block[0]->spin);
    //normalize(&refDirection);
    //printf("-------------------\n");
    //printf("the seed Orb is %d\n",block[0]->id);
    //printf("trial normal direction %.3f %.3f %.3f\n",refDirection.coor[0],refDirection.coor[1],refDirection.coor[2]);
    //double effectiveJ=diagonalDot(block[0]->linkStrength[0], refDirection, refDirection);
    //printf("the 0 linking has strength %.3f %.3f %.3f (original) %.3f (effective)\n",block[0]->linkStrength[0].coor[0],block[0]->linkStrength[1].coor[1],block[0]->linkStrength[2].coor[2],effectiveJ);
    while (expandBlock(p_beginIndex, p_endIndex, buffer, p_blockLen, block, *refDirection)==1)
    {
        // no code here
    }
    
    //printf("    Block size is %d\n",*p_blockLen);
    double tot_d_onsiteEnergy=0;
    for(i=0;i<*p_blockLen;i++){
        block[i]->isProjected=0;
        for (j = 0; j < block[i]->nlink; j++)
        {
            // exchange anisotropy
            block[i]->linkedOrb[j]->isProjected=0;
            Vec originalSj, Sj_parallel;
            equal(&originalSj, block[i]->linkedOrb[j]->spin);
            equal(&Sj_parallel,block[i]->linkedOrb[j]->transSpin);
            cTimes(&Sj_parallel,0.5);
            plusEqual(&originalSj,Sj_parallel);
            tot_d_onsiteEnergy+=block[i]->sDotN*diagonalDot(originalSj,block[i]->linkStrength[j],*refDirection);
        }
        // single-ion anisotropy
        tot_d_onsiteEnergy+=getDeltaOnsiteEnergy(block[i]);
        
    }
    for(i=0;i<*p_blockLen;i++) block[i]->inBlock=0;
    free(refDirection);
    // process the onsite anisotropy
    //printf("tot_d_onsiteEnergy=%.6f\n",tot_d_onsiteEnergy);
    if(tot_d_onsiteEnergy<=0 || exp(-tot_d_onsiteEnergy)>rand()/(double) RAND_MAX){
        for(i=0;i<*p_blockLen;i++){
            plusEqual(&block[i]->spin, block[i]->transSpin);
            //printf("    after update orb %d spin converted to %.3f %.3f %.3f\n",block[i]->id,block[i]->spin.coor[0],block[i]->spin.coor[1],block[i]->spin.coor[2]);
            //block[i]->inBlock=0;
            plusEqual(p_totSpin,block[i]->transSpin);
        }
        // update energy
        *p_energy=0.;
        for(i=0;i<totOrbs;i++) *p_energy+=getCorrEnergy(lattice+i);
        *p_energy/=2; // bond term
        for(i=0;i<totOrbs;i++) *p_energy+=getOnsiteEnergy(lattice+i); // onsite term
    }
}

void localUpdate(int totOrbs, Orb lattice[], double *p_energy, Vec *p_totSpin){
    //printf("start local updating\n");
    int seedID=rand()%totOrbs;  // chose one orb
    //printf("considering %d orb, its spin is %.3f %.3f %.3f\n",
    //         seedID,lattice[seedID].spin.coor[0],lattice[seedID].spin.coor[1],lattice[seedID].spin.coor[2]);
    Vec *refDirection=generateRandomVec(); // chose new direction
    //printf("try new spin direction, ref: %.3f %.3f %.3f\n",
    //       refDirection.coor[0],refDirection.coor[1],refDirection.coor[2]);
    double s1n=-2*dot(lattice[seedID].spin,*refDirection);
    //printf("projection s1n: %.3f\n",s1n);
    equal(&lattice[seedID].transSpin,*refDirection);
    cTimes(&lattice[seedID].transSpin,s1n);
    double corr=getDeltaCorrEnergy(lattice+seedID);
    corr+=getDeltaOnsiteEnergy(lattice+seedID);
    
    //printf("lead to the translation spin vector: %.3f %.3f %.3f and delta Ecorr: %.3f\n",
    //      lattice[seedID].transSpin.coor[0],lattice[seedID].transSpin.coor[1],lattice[seedID].transSpin.coor[2],corr);
    
    if(corr<=0 || exp(-corr)>rand()/(double) RAND_MAX){  // new direction is energertically favoured thus accept directly
        plusEqual(p_totSpin,lattice[seedID].transSpin);
        plusEqual(&lattice[seedID].spin,lattice[seedID].transSpin);
        *p_energy+=corr;
        //printf("since new direction is energertically lowerd thus we accept, energy: %.3f\n",*p_energy);
    }
    free(refDirection);
    return;
}

PyObject * blockUpdateMC(int totOrbs, double initSpin[totOrbs], double initD[totOrbs][3], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking][3], int linkedOrb[totOrbs][maxNLinking],
                   int ninterval, int nLat, int corrOrbPair[nLat][2], double flunc, double h,
                   int totOrb_rnorm, int nOrbInCluster, int rOrb[totOrb_rnorm], int rOrbCluster[totOrb_rnorm][nOrbInCluster], int linkedOrb_rnorm[totOrb_rnorm][maxNLinking]){
    // initialize lattice
    Orb lattice[totOrbs];
    //printf("hello here is C lib\n");
    establishLattice(lattice, totOrbs, initSpin, initD, h, flunc, maxNLinking, nlink, linkStrength, totOrb_rnorm, nOrbInCluster, rOrb, rOrbCluster);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb, totOrb_rnorm, rOrb, linkedOrb_rnorm);

    // initialize measurement
    double energy=0;
    double *p_energy=&energy;
    for(int i=0;i<totOrbs;i++)*p_energy+=getCorrEnergy(lattice+i);
    *p_energy/=2; // double counting
    for(int i=0;i<totOrbs;i++)*p_energy+=getOnsiteEnergy(lattice+i);
    //printf("gournd energy=%.3f nLat=%d\n",*p_energy, nLat);

    Vec totSpin;
    totSpin.x=0;totSpin.y=0;totSpin.z=0;
    Vec*p_totSpin=&totSpin;
    for(int i=0;i<totOrbs;i++) plusEqual(p_totSpin, lattice[i].spin);
    
    // initialize block
    for(int i=0;i<totOrbs;i++) {lattice[i].inBlock=0;lattice[i].isProjected=0;};

    for(int i=0;i<nthermal*ninterval;i++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin); //thermalization

    // printf("start sweeping\n");
    Vec spin_i, spin_i_r;
    Vec spin_j, spin_j_r;
    spin_i.x=0;spin_i.y=0;spin_i.z=0;spin_i_r.x=0;spin_i_r.y=0;spin_i_r.z=0;
    spin_j.x=0;spin_j.y=0;spin_j.z=0;spin_j_r.x=0;spin_j_r.y=0;spin_j_r.z=0;
    double spin_ij, spin_ij_r;
    spin_ij=0;spin_ij_r=0;

    double totEnergy=0, totEnergy_r=0;
    double E2=0, E2_r=0;
    double M=0,M2=0,M4=0;
    double M_tmp=0,MdotM_tmp=0,M_tot=0;

    Vec spin_direction;
    double spin_i_z, spin_j_z, spin_tot_z; // spin projected to net-spin
    double spin_i_h, spin_j_h, spin_tot_h; // spin projected to z aixis
    spin_i_z=0;spin_j_z=0;spin_tot_z=0;
    spin_i_h=0;spin_j_h=0;spin_tot_h=0;

    for(int i=0;i<nsweep;i++){
        for(int j=0;j<ninterval;j++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin);
        // find the main axis
        spin_direction.x=p_totSpin->x;
        spin_direction.y=p_totSpin->y;
        spin_direction.z=p_totSpin->z;
        normalize(&spin_direction);

        // spin statistics over space in each frame
        Vec spin_i_avg;
        Vec spin_j_avg;
        spin_i_avg.x=0;
        spin_i_avg.y=0;
        spin_i_avg.z=0;
        spin_j_avg.x=0;
        spin_j_avg.y=0;
        spin_j_avg.z=0;
        double spin_ij_avg=0.0;

        double spin_i_z_avg, spin_j_z_avg;
        double spin_i_h_avg, spin_j_h_avg;
        spin_i_z_avg=0;spin_j_z_avg=0;
        spin_i_h_avg=0;spin_j_h_avg=0;
        for(int j=0;j<nLat;j++){
            plusEqual(&spin_i_avg, lattice[corrOrbPair[j][0]].spin);
            plusEqual(&spin_j_avg, lattice[corrOrbPair[j][1]].spin);
            spin_ij_avg+=dot(lattice[corrOrbPair[j][0]].spin,lattice[corrOrbPair[j][1]].spin);
            
            // spin along main axis
            spin_i_z_avg+=dot(spin_direction,lattice[corrOrbPair[j][0]].spin);
            spin_j_z_avg+=dot(spin_direction,lattice[corrOrbPair[j][1]].spin);

            // spin projected to z axis
            spin_i_h_avg+=lattice[corrOrbPair[j][0]].spin.z;
            spin_j_h_avg+=lattice[corrOrbPair[j][1]].spin.z;
        }

        spin_i_z+=spin_i_z_avg/nLat;
        spin_j_z+=spin_j_z_avg/nLat;
        spin_tot_z+=(dot(spin_direction,*p_totSpin)/nLat);
        if(h<0.00001){// avoid faults time reversal symmetry
            spin_i_h+=fabs(spin_i_h_avg)/nLat;
            spin_j_h+=fabs(spin_j_h_avg)/nLat;
            spin_tot_h+=fabs(p_totSpin->z/nLat);
        }else{
            spin_i_h+=spin_i_h_avg/nLat;
            spin_j_h+=spin_j_h_avg/nLat;
            spin_tot_h+=p_totSpin->z/nLat;
        }


        M=sqrt(dot(*p_totSpin,*p_totSpin))/nLat;
        M2+=M*M;
        M4+=M*M*M*M;
        //calc auto-correlation
        M_tot+=M;
        MdotM_tmp+=M_tmp*M;
        M_tmp=M;

        cDivides(&spin_i_avg, nLat);
        cDivides(&spin_j_avg, nLat);
        vabs(&spin_i_avg);
        vabs(&spin_j_avg);
        plusEqual(&spin_i,spin_i_avg);
        plusEqual(&spin_j,spin_j_avg);
        spin_ij+=spin_ij_avg/nLat;

        double e_avg=*p_energy/totOrbs;
        totEnergy+=e_avg;
        E2+=e_avg*e_avg;

        // *************** statistics on renormalized lattice
        // spin statistics over space in each frame
        Vec spin_i_r_avg, spin_j_r_avg;
        spin_i_r_avg.x=0;spin_i_r_avg.y=0;spin_i_r_avg.z=0;
        spin_j_r_avg.x=0;spin_j_r_avg.y=0;spin_j_r_avg.z=0;
        int chosed_spin_i=0;
        int chosed_spin_j=0;
        int chosed_spin_ij=0;
        double spin_ij_r_avg=0.0;
        for(int j=0;j<nLat;j++){
            Orb *_orb_i=lattice+corrOrbPair[j][0];
            Orb *_orb_j=lattice+corrOrbPair[j][1];
            if(_orb_i->chosen>0){
                chosed_spin_i+=1;
                //printf("_orb_i.spin: %.3f %.3f\n",_orb_i->spin.x,_orb_i->spin.y);
                //Vec majority=getMajoritySpin(_orb_i);
                //printf("majority: %.3f %.3f\n",majority.x,majority.y);
                plusEqual(&spin_i_r_avg, getMajoritySpin(_orb_i));
            }
            if(_orb_j->chosen>0){
                //printf("_orb_j.spin: %.3f %.3f\n",_orb_j->spin.x,_orb_j->spin.y);
                //Vec majority=getMajoritySpin(_orb_j);
                //printf("majority: %.3f %.3f\n",majority.x,majority.y);
                chosed_spin_j+=1;
                plusEqual(&spin_j_r_avg, getMajoritySpin(_orb_j));
            }
            if(_orb_i->chosen*_orb_j->chosen>0){
                chosed_spin_ij+=1;
                //Vec mi=getMajoritySpin(_orb_i);
                //Vec mj=getMajoritySpin(_orb_j);
                //printf("dot product: %.3f,%.3f dot %.3f,%.3f = %.3f\n",mi.x,mi.y,mj.x,mj.y,mi.x*mj.x+mi.y*mj.y);
                spin_ij_r_avg+=dot(getMajoritySpin(_orb_i),getMajoritySpin(_orb_j));
            }
            
        }

        //printf("chosen orb quantity: %d\n",chosed_spin_i);
        //printf("accumulated spin %.3f %.3f\n",spin_i_r_avg.x,spin_i_r_avg.y);
        cDivides(&spin_i_r_avg, chosed_spin_i);
        cDivides(&spin_j_r_avg, chosed_spin_j);
        vabs(&spin_i_r_avg);
        vabs(&spin_j_r_avg);
        plusEqual(&spin_i_r,spin_i_r_avg);
        plusEqual(&spin_j_r,spin_j_r_avg);
        spin_ij_r+=spin_ij_r_avg/chosed_spin_ij;

        // energy statistics
        double e_avg_rnorm=0;
        for(int j=0;j<totOrbs;j++){ // calc. bond energy in renormalized system
            if(lattice[j].chosen>0) e_avg_rnorm+=getCorrEnergy_rnorm(lattice+j);
        }
        e_avg_rnorm/=2; // double counting
        for(int j=0;j<totOrbs;j++){ // onsite part
            if(lattice[j].chosen>0) e_avg_rnorm+=getOnsiteEnergy(lattice+j);
        }
        e_avg_rnorm/=totOrb_rnorm;
        totEnergy_r+=e_avg_rnorm;
        E2_r+=e_avg_rnorm*e_avg_rnorm;
    }
    double U4=(M2/nsweep)*(M2/nsweep)/(M4/nsweep);
    double autoCorr=(MdotM_tmp/nsweep-M_tot/nsweep*M_tot/nsweep);
    PyObject *Data;
    Data=PyTuple_New(26);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i.x/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_i.y/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(spin_i.z/nsweep));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(spin_j.x/nsweep));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(spin_j.y/nsweep));
    PyTuple_SetItem(Data, 5, PyFloat_FromDouble(spin_j.z/nsweep));
    PyTuple_SetItem(Data, 6, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 7, PyFloat_FromDouble(autoCorr));
    PyTuple_SetItem(Data, 8, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 9, PyFloat_FromDouble(E2/nsweep));
    PyTuple_SetItem(Data,10, PyFloat_FromDouble(U4));
    PyTuple_SetItem(Data,11, PyFloat_FromDouble(spin_i_r.x/nsweep));
    PyTuple_SetItem(Data,12, PyFloat_FromDouble(spin_i_r.y/nsweep));
    PyTuple_SetItem(Data,13, PyFloat_FromDouble(spin_i_r.z/nsweep));
    PyTuple_SetItem(Data,14, PyFloat_FromDouble(spin_j_r.x/nsweep));
    PyTuple_SetItem(Data,15, PyFloat_FromDouble(spin_j_r.y/nsweep));
    PyTuple_SetItem(Data,16, PyFloat_FromDouble(spin_j_r.z/nsweep));
    PyTuple_SetItem(Data,17, PyFloat_FromDouble(spin_ij_r/nsweep));
    PyTuple_SetItem(Data,18, PyFloat_FromDouble(totEnergy_r/nsweep));
    PyTuple_SetItem(Data,19, PyFloat_FromDouble(E2_r/nsweep));
    PyTuple_SetItem(Data,20, PyFloat_FromDouble(spin_i_z/nsweep));
    PyTuple_SetItem(Data,21, PyFloat_FromDouble(spin_j_z/nsweep));
    PyTuple_SetItem(Data,22, PyFloat_FromDouble(spin_tot_z/nsweep));
    PyTuple_SetItem(Data,23, PyFloat_FromDouble(spin_i_h/nsweep));
    PyTuple_SetItem(Data,24, PyFloat_FromDouble(spin_j_h/nsweep));
    PyTuple_SetItem(Data,25, PyFloat_FromDouble(spin_tot_h/nsweep));
    return Data;
}


// self.totOrbs, initSpin, nthermal, nsweep, maxNLinking, nlinking, linkStrength, linkData
PyObject * localUpdateMC(int totOrbs, double initSpin[totOrbs], double initD[totOrbs][3], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking][3], int linkedOrb[totOrbs][maxNLinking],
                   int ninterval, int nLat, int corrOrbPair[nLat][2], double flunc, double h,
                   int totOrb_rnorm, int nOrbInCluster, int rOrb[totOrb_rnorm], int rOrbCluster[totOrb_rnorm][nOrbInCluster], int linkedOrb_rnorm[totOrb_rnorm][maxNLinking]){
    // initialize lattice
    Orb lattice[totOrbs];
    //printf("hello here is C lib\n");
    establishLattice(lattice, totOrbs, initSpin, initD, h, flunc, maxNLinking, nlink, linkStrength, totOrb_rnorm, nOrbInCluster, rOrb, rOrbCluster);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb, totOrb_rnorm, rOrb, linkedOrb_rnorm);

    // initialize measurement
    double energy=0;
    double *p_energy=&energy;
    for(int i=0;i<totOrbs;i++)*p_energy+=getCorrEnergy(lattice+i);
    *p_energy/=2; // double counting for bond term
    for(int i=0;i<totOrbs;i++)*p_energy+=getOnsiteEnergy(lattice+i); // onsite term

    Vec totSpin;
    totSpin.x=0;totSpin.y=0;totSpin.z=0;
    Vec*p_totSpin=&totSpin;
    for(int i=0;i<totOrbs;i++) {
        //printf("%.3f %.3f %.3f\n",lattice[i].spin.x,lattice[i].spin.y,lattice[i].spin.y);
        plusEqual(&totSpin, lattice[i].spin);
    }
    
    //localUpdate(totOrbs, lattice, p_energy, p_totSpin);
    //printf("initial total spin: %.3f %.3f %.3f energy: %.3f\n",totSpin.x,totSpin.y,totSpin.z,*p_energy);
    for(int i=0;i<ninterval*nthermal;i++) localUpdate(totOrbs, lattice, p_energy, p_totSpin); //thermalization

    Vec spin_i, spin_i_r;
    Vec spin_j, spin_j_r;
    spin_i.x=0;spin_i.y=0;spin_i.z=0;spin_i_r.x=0;spin_i_r.y=0;spin_i_r.z=0;
    spin_j.x=0;spin_j.y=0;spin_j.z=0;spin_j_r.x=0;spin_j_r.y=0;spin_j_r.z=0;
    double spin_ij, spin_ij_r;
    spin_ij=0;spin_ij_r=0;

    double totEnergy=0, totEnergy_r=0;
    double E2=0, E2_r=0;
    double M=0,M2=0,M4=0;
    double M_tmp=0,MdotM_tmp=0,M_tot=0;

    Vec spin_direction;
    double spin_i_z, spin_j_z, spin_tot_z;
    double spin_i_h, spin_j_h, spin_tot_h;
    spin_i_z=0;spin_j_z=0;spin_tot_z=0;
    spin_i_h=0;spin_j_h=0;spin_tot_h=0;

    for(int i=0;i<nsweep;i++){
        for(int j=0;j<ninterval;j++) localUpdate(totOrbs, lattice, p_energy, p_totSpin);
        // find the main axis
        spin_direction.x=p_totSpin->x;
        spin_direction.y=p_totSpin->y;
        spin_direction.z=p_totSpin->z;
        normalize(&spin_direction);
        

        // spin statistics over space in each frame
        Vec spin_i_avg;
        Vec spin_j_avg;
        spin_i_avg.x=0;
        spin_i_avg.y=0;
        spin_i_avg.z=0;
        spin_j_avg.x=0;
        spin_j_avg.y=0;
        spin_j_avg.z=0;
        double spin_ij_avg=0.0;

        double spin_i_z_avg, spin_j_z_avg;
        double spin_i_h_avg, spin_j_h_avg;
        spin_i_z_avg=0;spin_j_z_avg=0;
        spin_i_h_avg=0;spin_j_h_avg=0;
        for(int j=0;j<nLat;j++){
            plusEqual(&spin_i_avg, lattice[corrOrbPair[j][0]].spin);
            plusEqual(&spin_j_avg, lattice[corrOrbPair[j][1]].spin);
            spin_ij_avg+=dot(lattice[corrOrbPair[j][0]].spin,lattice[corrOrbPair[j][1]].spin);

            // spin along main axis
            spin_i_z_avg+=dot(spin_direction,lattice[corrOrbPair[j][0]].spin);
            spin_j_z_avg+=dot(spin_direction,lattice[corrOrbPair[j][1]].spin);

            // spin projected to z axis
            spin_i_h_avg+=lattice[corrOrbPair[j][0]].spin.z;
            spin_j_h_avg+=lattice[corrOrbPair[j][1]].spin.z;
        }

        spin_i_z+=spin_i_z_avg/nLat;
        spin_j_z+=spin_j_z_avg/nLat;
        spin_tot_z+=(dot(spin_direction,*p_totSpin)/nLat);
        if(h<0.00001){// avoid faults time reversal symmetry
            spin_i_h+=fabs(spin_i_h_avg)/nLat;
            spin_j_h+=fabs(spin_j_h_avg)/nLat;
            spin_tot_h+=fabs(p_totSpin->z)/nLat;
        }else{
            spin_i_h+=spin_i_h_avg/nLat;
            spin_j_h+=spin_j_h_avg/nLat;
            spin_tot_h+=p_totSpin->z/nLat;
        }

        M=sqrt(dot(spin_i_avg,spin_i_avg))/nLat;
        M2+=M*M;
        M4+=M*M*M*M;
        //calc auto-correlation
        M_tot+=M;
        MdotM_tmp+=M_tmp*M;
        M_tmp=M;

        cDivides(&spin_i_avg, nLat);
        cDivides(&spin_j_avg, nLat);
        vabs(&spin_i_avg);
        vabs(&spin_j_avg);
        plusEqual(&spin_i,spin_i_avg);
        plusEqual(&spin_j,spin_j_avg);
        spin_ij+=spin_ij_avg/nLat;


        // save energy variables
        double e_avg=*p_energy/totOrbs;
        totEnergy+=e_avg;
        E2+=e_avg*e_avg;

        // *************** statistics on renormalized lattice
        // spin statistics over space in each frame
        Vec spin_i_r_avg, spin_j_r_avg;
        spin_i_r_avg.x=0;spin_i_r_avg.y=0;spin_i_r_avg.z=0;
        spin_j_r_avg.x=0;spin_j_r_avg.y=0;spin_j_r_avg.z=0;
        int chosed_spin_i=0;
        int chosed_spin_j=0;
        int chosed_spin_ij=0;
        double spin_ij_r_avg=0.0;
        for(int j=0;j<nLat;j++){
            Orb *_orb_i=lattice+corrOrbPair[j][0];
            Orb *_orb_j=lattice+corrOrbPair[j][1];
            if(_orb_i->chosen>0){
                chosed_spin_i+=1;
                //printf("_orb_i.spin: %.3f %.3f\n",_orb_i->spin.x,_orb_i->spin.y);
                //Vec majority=getMajoritySpin(_orb_i);
                //printf("majority: %.3f %.3f\n",majority.x,majority.y);
                plusEqual(&spin_i_r_avg, getMajoritySpin(_orb_i));
            }
            if(_orb_j->chosen>0){
                //printf("_orb_j.spin: %.3f %.3f\n",_orb_j->spin.x,_orb_j->spin.y);
                //Vec majority=getMajoritySpin(_orb_j);
                //printf("majority: %.3f %.3f\n",majority.x,majority.y);
                chosed_spin_j+=1;
                plusEqual(&spin_j_r_avg, getMajoritySpin(_orb_j));
            }
            if(_orb_i->chosen*_orb_j->chosen>0){
                chosed_spin_ij+=1;
                //Vec mi=getMajoritySpin(_orb_i);
                //Vec mj=getMajoritySpin(_orb_j);
                //printf("dot product: %.3f,%.3f dot %.3f,%.3f = %.3f\n",mi.x,mi.y,mj.x,mj.y,mi.x*mj.x+mi.y*mj.y);
                spin_ij_r_avg+=dot(getMajoritySpin(_orb_i),getMajoritySpin(_orb_j));
            }
            
        }

        //printf("chosen orb quantity: %d\n",chosed_spin_i);
        //printf("accumulated spin %.3f %.3f\n",spin_i_r_avg.x,spin_i_r_avg.y);
        cDivides(&spin_i_r_avg, chosed_spin_i);
        cDivides(&spin_j_r_avg, chosed_spin_j);
        vabs(&spin_i_r_avg);
        vabs(&spin_j_r_avg);
        plusEqual(&spin_i_r,spin_i_r_avg);
        plusEqual(&spin_j_r,spin_j_r_avg);
        spin_ij_r+=spin_ij_r_avg/chosed_spin_ij;

        // energy statistics
        double e_avg_rnorm=0;
        for(int j=0;j<totOrbs;j++){ // calc. bond energy in renormalized system
            if(lattice[j].chosen>0) e_avg_rnorm+=getCorrEnergy_rnorm(lattice+j);
        }
        e_avg_rnorm/=2; // double counting
        for(int j=0;j<totOrbs;j++){ // onsite part
            if(lattice[j].chosen>0) e_avg_rnorm+=getOnsiteEnergy(lattice+j);
        }
        e_avg_rnorm/=totOrb_rnorm;
        totEnergy_r+=e_avg_rnorm;
        E2_r+=e_avg_rnorm*e_avg_rnorm;
    }
    double U4=(M2/nsweep)*(M2/nsweep)/(M4/nsweep);
    double autoCorr=(MdotM_tmp/nsweep-M_tot/nsweep*M_tot/nsweep);
    PyObject *Data;
    Data=PyTuple_New(26);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i.x/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_i.y/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(spin_i.z/nsweep));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(spin_j.x/nsweep));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(spin_j.y/nsweep));
    PyTuple_SetItem(Data, 5, PyFloat_FromDouble(spin_j.z/nsweep));
    PyTuple_SetItem(Data, 6, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 7, PyFloat_FromDouble(autoCorr));
    PyTuple_SetItem(Data, 8, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 9, PyFloat_FromDouble(E2/nsweep));
    PyTuple_SetItem(Data,10, PyFloat_FromDouble(U4));
    PyTuple_SetItem(Data,11, PyFloat_FromDouble(spin_i_r.x/nsweep));
    PyTuple_SetItem(Data,12, PyFloat_FromDouble(spin_i_r.y/nsweep));
    PyTuple_SetItem(Data,13, PyFloat_FromDouble(spin_i_r.z/nsweep));
    PyTuple_SetItem(Data,14, PyFloat_FromDouble(spin_j_r.x/nsweep));
    PyTuple_SetItem(Data,15, PyFloat_FromDouble(spin_j_r.y/nsweep));
    PyTuple_SetItem(Data,16, PyFloat_FromDouble(spin_j_r.z/nsweep));
    PyTuple_SetItem(Data,17, PyFloat_FromDouble(spin_ij_r/nsweep));
    PyTuple_SetItem(Data,18, PyFloat_FromDouble(totEnergy_r/nsweep));
    PyTuple_SetItem(Data,19, PyFloat_FromDouble(E2_r/nsweep));
    PyTuple_SetItem(Data,20, PyFloat_FromDouble(spin_i_z/nsweep));
    PyTuple_SetItem(Data,21, PyFloat_FromDouble(spin_j_z/nsweep));
    PyTuple_SetItem(Data,22, PyFloat_FromDouble(spin_tot_z/nsweep));
    PyTuple_SetItem(Data,23, PyFloat_FromDouble(spin_i_h/nsweep));
    PyTuple_SetItem(Data,24, PyFloat_FromDouble(spin_j_h/nsweep));
    PyTuple_SetItem(Data,25, PyFloat_FromDouble(spin_tot_h/nsweep));
    return Data;
}