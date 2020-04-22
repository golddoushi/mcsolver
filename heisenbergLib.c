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
Vec *generateRandomVec(){
    Vec *direction=(Vec*)malloc(sizeof(Vec));
    direction->x=rand()/32767.0-0.5;
    direction->y=rand()/32767.0-0.5;
    direction->z=rand()/32767.0-0.5;
    double len2=dot(*direction, *direction);
    if (len2>0.5)
    {
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
    Vec spin;
    Vec transSpin;
    int nlink;
    Vec *linkStrength;
    int inBlock;
    struct Orb **linkedOrb;
    Vec onsiteAnisotropy;

    double d_onsiteEnergy;
    double sDotN;
    int isProjected;
    
}Orb;

//establishLattice(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength);
void establishLattice(Orb *lattice, int totOrbs, double initSpin[totOrbs], double initD[totOrbs][3], double flunc, int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking][3]){
    //printf("establishing whole lattice with %d orbs and %d linkings for each orb\n",totOrbs,maxNLinking);
    for(int i=0;i<totOrbs;i++){
        //printf("check point-1, entering setting %d",i);
        lattice[i].id=i;
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
    }
}

void establishLinking(Orb *lattice, int totOrbs, int maxNLinking, int nlink[totOrbs], int linkedOrb[totOrbs][maxNLinking]){
    for(int iorb=0;iorb<totOrbs;iorb++){
        lattice[iorb].linkedOrb=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb][ilink];
        }
    }
}

double getCorrEnergy(Orb *source){
    double corr=0;
    for(int i=0;i<source->nlink;i++){
        corr+=diagonalDot(source->linkStrength[i],source->spin,source->linkedOrb[i]->spin);
    }
    return corr;
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
    return source->onsiteAnisotropy.x*(s1x*s1x-source->spin.x*source->spin.x)+
           source->onsiteAnisotropy.y*(s1y*s1y-source->spin.y*source->spin.y)+
           source->onsiteAnisotropy.z*(s1z*s1z-source->spin.z*source->spin.z);
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
        outOrb->isProjected==1;
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
            if(corr<0 && (1-exp(corr))>rand()/32767.0){
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
    if(tot_d_onsiteEnergy<=0 || exp(-tot_d_onsiteEnergy)>rand()/32767.0){
        for(i=0;i<*p_blockLen;i++){
            plusEqual(&block[i]->spin, block[i]->transSpin);
            //printf("    after update orb %d spin converted to %.3f %.3f %.3f\n",block[i]->id,block[i]->spin.coor[0],block[i]->spin.coor[1],block[i]->spin.coor[2]);
            //block[i]->inBlock=0;
            plusEqual(p_totSpin,block[i]->transSpin);
        }
        *p_energy=0.;
        for(i=0;i<totOrbs;i++){
            *p_energy+=getCorrEnergy(lattice+i);
        }
        *p_energy/=2;
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
    
    if(corr<=0 || exp(-corr)>rand()/32767.0){  // new direction is energertically favoured thus accept directly
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
                   int ninterval, int nLat, int corrOrbPair[nLat][2], double flunc, double h){
    // initialize lattice
    Orb lattice[totOrbs];
    //printf("hello here is C lib\n");
    establishLattice(lattice, totOrbs, initSpin, initD, flunc, maxNLinking, nlink, linkStrength);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb);

    // initialize measurement
    double energy=0;
    double *p_energy=&energy;
    Vec totSpin;
    totSpin.x=0;totSpin.y=0;totSpin.z=0;
    Vec*p_totSpin=&totSpin;
    for(int i=0;i<totOrbs;i++) plusEqual(p_totSpin, lattice[i].spin);
    
    // initialize block
    for(int i=0;i<totOrbs;i++) {lattice[i].inBlock=0;lattice[i].isProjected=0;};

    for(int i=0;i<nthermal*ninterval;i++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin); //thermalization

    // printf("start sweeping\n");
    Vec spin_i;
    Vec spin_j;
    spin_i.x=0;
    spin_i.y=0;
    spin_i.z=0;
    spin_j.x=0;
    spin_j.y=0;
    spin_j.z=0;
    double spin_ij;
    double totEnergy=0;
    double E2=0;

    for(int i=0;i<nsweep;i++){
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
        for(int j=0;j<nLat;j++){
            plusEqual(&spin_i_avg, lattice[corrOrbPair[j][0]].spin);
            plusEqual(&spin_j_avg, lattice[corrOrbPair[j][1]].spin);
            spin_ij_avg+=dot(lattice[corrOrbPair[j][0]].spin,lattice[corrOrbPair[j][1]].spin);
            
        }
        cDivides(&spin_i_avg, nLat);
        cDivides(&spin_j_avg, nLat);
        vabs(&spin_i_avg);
        vabs(&spin_j_avg);
        plusEqual(&spin_i,spin_i_avg);
        plusEqual(&spin_j,spin_j_avg);
        spin_ij+=spin_ij_avg/nLat;

        double e_avg=*p_energy/nLat;
        totEnergy+=e_avg;
        E2+=e_avg*e_avg;

        for(int j=0;j<ninterval;j++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin);
    }
    PyObject *Data;
    Data=PyTuple_New(9);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i.x/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_i.y/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(spin_i.z/nsweep));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(spin_j.x/nsweep));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(spin_j.y/nsweep));
    PyTuple_SetItem(Data, 5, PyFloat_FromDouble(spin_j.z/nsweep));
    PyTuple_SetItem(Data, 6, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 7, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 8, PyFloat_FromDouble(E2/nsweep));
    return Data;
}


// self.totOrbs, initSpin, nthermal, nsweep, maxNLinking, nlinking, linkStrength, linkData
PyObject * localUpdateMC(int totOrbs, double initSpin[totOrbs], double initD[totOrbs][3], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking][3], int linkedOrb[totOrbs][maxNLinking],
                   int ninterval, int nLat, int corrOrbPair[nLat][2], double flunc, double h){
    // initialize lattice
    Orb lattice[totOrbs];
    //printf("hello here is C lib\n");
    establishLattice(lattice, totOrbs, initSpin, initD, flunc, maxNLinking, nlink, linkStrength);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb);

    // initialize measurement
    double energy=0;
    double *p_energy=&energy;
    for(int i=0;i<totOrbs;i++){
        *p_energy+=getCorrEnergy(lattice+i);
    }
    *p_energy/=2;
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

    Vec spin_i;
    Vec spin_j;
    spin_i.x=0;
    spin_i.y=0;
    spin_i.z=0;
    spin_j.x=0;
    spin_j.y=0;
    spin_j.z=0;
    double spin_ij;
    double totEnergy=0;
    double E2=0;

    for(int i=0;i<nsweep;i++){
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
        for(int j=0;j<nLat;j++){
            plusEqual(&spin_i_avg, lattice[corrOrbPair[j][0]].spin);
            plusEqual(&spin_j_avg, lattice[corrOrbPair[j][1]].spin);
            spin_ij_avg+=dot(lattice[corrOrbPair[j][0]].spin,lattice[corrOrbPair[j][1]].spin);
            
        }
        cDivides(&spin_i_avg, nLat);
        cDivides(&spin_j_avg, nLat);
        vabs(&spin_i_avg);
        vabs(&spin_j_avg);
        plusEqual(&spin_i,spin_i_avg);
        plusEqual(&spin_j,spin_j_avg);
        spin_ij+=spin_ij_avg/nLat;

        double e_avg=*p_energy/nLat;
        totEnergy+=e_avg;
        E2+=e_avg*e_avg;

        for(int j=0;j<ninterval;j++) localUpdate(totOrbs, lattice, p_energy, p_totSpin);
    }
    
    PyObject *Data;
    Data=PyTuple_New(9);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i.x/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_i.y/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(spin_i.z/nsweep));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(spin_j.x/nsweep));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(spin_j.y/nsweep));
    PyTuple_SetItem(Data, 5, PyFloat_FromDouble(spin_j.z/nsweep));
    PyTuple_SetItem(Data, 6, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 7, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 8, PyFloat_FromDouble(E2/nsweep));
    return Data;
}