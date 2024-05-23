#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Vec
{
    int dimension;
    double x,y;
}Vec;
typedef struct Vec4
{
    int dimension;
    double xx,yy,xy,yx;
}Vec4;
// vec1=abs(vec1)
void normalize(Vec *vec1){
    double len=sqrt(vec1->x*vec1->x+vec1->y*vec1->y);
    if (len<1e-5) return;
    vec1->x/=len;
    vec1->y/=len;
}
// vec1=vec2
void equal(Vec *vec1, Vec vec2){
    vec1->x=vec2.x;
    vec1->y=vec2.y;
}
// vec1*=c
void cTimes(Vec *vec1, double c){
    vec1->x*=c;
    vec1->y*=c;
}
// vec1=fabs(vec1)
void vabs(Vec *vec1){
    vec1->x=fabs(vec1->x);
    vec1->y=fabs(vec1->y);
}
// vec1/=c
void cDivides(Vec *vec1, double c){
    vec1->x/=c;
    vec1->y/=c;
}
// vec1+=vec2
void plusEqual(Vec *vec1, Vec vec2){
    vec1->x+=vec2.x;
    vec1->y+=vec2.y;
}
// vec1-=vec2
void minusEqual(Vec *vec1, Vec vec2){
    vec1->x-=vec2.x;
    vec1->y-=vec2.y;
}
// vec1 dot vec2
double dot(Vec vec1, Vec vec2){
    return vec1.x*vec2.x+vec1.y*vec2.y;
}

double diagonalDot(Vec vec1, Vec vec2, Vec4 matrix){
    return vec1.x*vec2.x*matrix.xx+
           vec1.y*vec2.y*matrix.yy+
           vec1.x*vec2.y*matrix.xy+
           vec1.y*vec2.x*matrix.yx;
}
double diagonalDot_simple(Vec vec1, Vec vec2, Vec4 matrix){
    return vec1.x*vec2.x*matrix.xx+
           vec1.y*vec2.y*matrix.yy;
}
double (*p_diagonalDot)(Vec vec1, Vec vec2, Vec4 matrix);

double gaussian_distr[RAND_MAX];

Vec *generateRandomVec(void){
    double x=rand()/(double) RAND_MAX-0.5;
    double y=rand()/(double) RAND_MAX-0.5;
    double len2=(x*x+y*y);
    if (len2>0.25)
    {
        return generateRandomVec();
    }else
    {
        double len=sqrt(len2);
        Vec *direction=(Vec*)malloc(sizeof(Vec));
        direction->x=x/len;
        direction->y=y/len;
        return direction;
    }
}

typedef struct Orb
{
    int id;
    double S; // maximum spin number
    Vec spin;
    Vec transSpin;
    Vec perpenSpin;
    int nlink;
    Vec4 *linkStrength;
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
void establishLattice(Orb *lattice, int totOrbs, double initSpin[totOrbs], double initD[totOrbs*3], double flunc, int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs*maxNLinking*9], double h,
                      int totOrb_rnorm, int nOrbInCluster, int rOrb[totOrb_rnorm], int rOrbCluster[totOrb_rnorm*nOrbInCluster]){
    //printf("establishing whole lattice with %d orbs and %d linkings for each orb\n",totOrbs,maxNLinking);
    int i;
    for(i=0;i<totOrbs;i++){
        //printf("check point-1, entering setting %d",i);
        lattice[i].id=i;
        lattice[i].S=initSpin[i];
        lattice[i].spin.x=initSpin[i];
        lattice[i].spin.y=0;

        Vec *fluncSpin=generateRandomVec();
        cTimes(fluncSpin, flunc);
        plusEqual(&lattice[i].spin,*fluncSpin);
        normalize(&lattice[i].spin);
        cTimes(&lattice[i].spin,fabs(initSpin[i]));
        free(fluncSpin);

        lattice[i].onsiteAnisotropy.x=initD[i*3+0];
        lattice[i].onsiteAnisotropy.y=initD[i*3+1];

        lattice[i].h=h;
        //printf("orb %d spin: %.3f %.3f %.3f\n",i,lattice[i].spin.coor[0],lattice[i].spin.coor[1],lattice[i].spin.coor[2]);
        lattice[i].transSpin.x=0;  // allocate trans spin vector for each orb
        lattice[i].transSpin.y=0;
        lattice[i].perpenSpin.x=0;
        lattice[i].perpenSpin.y=0;
        lattice[i].nlink=nlink[i];
        //printf("check point 1, orb: %d\n",lattice[i].id);
        lattice[i].linkStrength=(Vec4*)malloc(lattice[i].nlink*sizeof(Vec4)); // allocate strength for each linking
        //printf("check point 2, total links:%d\n",nlink[i]);
        for(int j=0;j<nlink[i];j++){
            //lattice[i].linkStrength[j].coor=(double*)malloc(3*sizeof(double));
            //printf("check point 3\n");
            lattice[i].linkStrength[j].xx=linkStrength[i*maxNLinking*9+j*9+0];
            lattice[i].linkStrength[j].yy=linkStrength[i*maxNLinking*9+j*9+1];
            lattice[i].linkStrength[j].xy=linkStrength[i*maxNLinking*9+j*9+3];
            lattice[i].linkStrength[j].yx=linkStrength[i*maxNLinking*9+j*9+6];
            //printf("check point 4, link:%d, strength: %.3f %.3f\n",j,lattice[i].linkStrength[j].x,lattice[i].linkStrength[j].y);
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
            lattice[id].orb_cluster[iorb]=lattice+rOrbCluster[i*nOrbInCluster+iorb];
            //printf("    from input id=%d orb%d\n",rOrbCluster[id][iorb],lattice[id].orb_cluster[iorb]->id);
        }
    }
    //printf("orbitals successfully built\n");
}

void establishLinking(Orb *lattice, int totOrbs, int maxNLinking, int nlink[totOrbs], int linkedOrb[totOrbs*maxNLinking],
                      int totOrb_rnorm, int rOrb[totOrb_rnorm], int linkedOrb_rnorm[totOrb_rnorm*maxNLinking]){
    for(int iorb=0;iorb<totOrbs;iorb++){
        lattice[iorb].linkedOrb=(Orb**)malloc(lattice[iorb].nlink*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb*maxNLinking+ilink];
        }
    }
    for(int i=0;i<totOrb_rnorm;i++){
        int iorb=rOrb[i];
        lattice[iorb].linkedOrb_rnorm=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb_rnorm[ilink]=lattice+linkedOrb_rnorm[i*maxNLinking+ilink];
        }
    }
    //printf("bonds successfully built\n");
}

double getCorrEnergy(Orb *source){
    double corr=0;
    for(int i=0;i<source->nlink;i++){
        //printf("getCorrEnergy\n");
        corr+=(*p_diagonalDot)(source->spin,source->linkedOrb[i]->spin,source->linkStrength[i]);
    }
    return corr;
}

double getOnsiteEnergy(Orb *source){
    return source->onsiteAnisotropy.x*source->spin.x*source->spin.x+
           source->onsiteAnisotropy.y*source->spin.y*source->spin.y-
           source->h*source->spin.x;
}

Vec getMajoritySpin(Orb*_orb){  // core algorithm for renormalization
    //return _orb->spin; // decimation algorithm
    Vec avgSpin;
    avgSpin.x=0;
    avgSpin.y=0;
    for(int iorb=0;iorb<_orb->nOrbInCluster;iorb++){
        plusEqual(&avgSpin,_orb->orb_cluster[iorb]->spin);
    }
    //cDivides(&avgSpin,_orb->nOrbInCluster);
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
        corr+=(*p_diagonalDot)(avgSpin_source,avgSpin_target,source->linkStrength[i]);
    }
    //printf("Ecorr=%.3f\n",corr);
    return corr;
}

double getOnsiteEnergy_rnorm(Orb *source){
    Vec avgSpin_source=getMajoritySpin(source);
    return source->onsiteAnisotropy.x*avgSpin_source.x*avgSpin_source.x+
           source->onsiteAnisotropy.y*avgSpin_source.y*avgSpin_source.y-
           source->h*avgSpin_source.x;
}

double getDeltaCorrEnergy(Orb *source){
    double corr=0;
    for(int i=0;i<source->nlink;i++){
        //printf("getDeltaCorrEnergy\n");
        corr+=(*p_diagonalDot)(source->transSpin,source->linkedOrb[i]->spin,source->linkStrength[i]);
    }
    return corr;
}

double getDeltaOnsiteEnergy(Orb *source){
    double s1x=source->spin.x+source->transSpin.x;
    double s1y=source->spin.y+source->transSpin.y;
    return source->onsiteAnisotropy.x*(s1x*s1x-source->spin.x*source->spin.x)+
           source->onsiteAnisotropy.y*(s1y*s1y-source->spin.y*source->spin.y)-
           source->h*source->transSpin.x;
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
        outOrb->sDotN=-dot(outOrb->spin,refDirection);
        equal(&outOrb->transSpin, refDirection);
        cTimes(&outOrb->transSpin, outOrb->sDotN);
        equal(&outOrb->perpenSpin,outOrb->spin);
        plusEqual(&outOrb->perpenSpin,outOrb->transSpin);
        outOrb->isProjected=1;
    }
    //printf("center orb: %d\n", outOrb->id);
    //printf("projection along axis: %.3f\n", -s1n/2);
    //printf("variation of spin: %.3f %.3f\n",outOrb->transSpin.x,outOrb->transSpin.y);
    //printf("it has %d neighbor orbs\n",outOrb->nlink);
    int i;
    for(i=0;i<outOrb->nlink;i++){
        //double effectiveJ=diagonalDot(outOrb->linkStrength[i], refDirection, refDirection);
        //printf("the %d linking has strength %.3f %.3f (original) %.3f (effective)\n",i,outOrb->linkStrength[i].x,outOrb->linkStrength[i].y,effectiveJ);
        Orb *linkedOrb=outOrb->linkedOrb[i];
        //printf("consider to add orb %d\n",linkedOrb->id);
        //printf("      considering the %d orb which is linking to %d orb, it is %d in block \n", linkedOrb->id, outOrb->id, linkedOrb->inBlock);
        if(linkedOrb->inBlock==0){
            //printf("projection along axis: %.3f\n", s2n);
            //double corr=-s1n*diagonalDot(refDirection,outOrb->linkStrength[i],linkedOrb->spin); // bond strength
            //printf("      spin of orb %d is %.3f %.3f %.3f and bond strength is %.3f\n",linkedOrb->id,linkedOrb->spin.coor[0],linkedOrb->spin.coor[1],linkedOrb->spin.coor[2],corr);
            if (linkedOrb->isProjected==0)
            {
                linkedOrb->sDotN=-dot(linkedOrb->spin,refDirection);
                equal(&linkedOrb->transSpin, refDirection);
                cTimes(&linkedOrb->transSpin, linkedOrb->sDotN);
                equal(&linkedOrb->perpenSpin, linkedOrb->spin);
                plusEqual(&linkedOrb->perpenSpin, linkedOrb->transSpin);
                linkedOrb->isProjected=1;
            }
            //printf("229\n");
            double corr=2*outOrb->sDotN*linkedOrb->sDotN*(*p_diagonalDot)(refDirection,refDirection,outOrb->linkStrength[i]);
            
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
    for(int i=0;i<totOrbs;i++){lattice[i].isProjected=0;lattice[i].inBlock=0;} // initialize all orb status
    Orb **block=(Orb**)malloc(totOrbs*sizeof(Orb*));
    Orb **buffer=(Orb**)malloc(totOrbs*sizeof(Orb*));
    int seedID=(rand()*RAND_MAX+rand())%totOrbs;
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
    //printf("trial normal direction %.3f %.3f\n",refDirection->x,refDirection->y);
    //double effectiveJ=diagonalDot(block[0]->linkStrength[0], *refDirection, *refDirection);
    //printf("the 0 linking has strength %.3f %.3f (original) %.3f (effective)\n",block[0]->linkStrength[0].x,block[0]->linkStrength[1].y,effectiveJ);
    while (expandBlock(p_beginIndex, p_endIndex, buffer, p_blockLen, block, *refDirection)==1)
    {
        // no code here
    }
    
    //printf("    Block size is %d\n",*p_blockLen);
    double tot_d_onsiteEnergy=0;
    // exchange anisotropy
    for(i=0;i<*p_blockLen;i++){
        for (j = 0; j < block[i]->nlink; j++){
            double source_anisotropy=block[i]->sDotN*(*p_diagonalDot)(*refDirection,block[i]->linkedOrb[j]->perpenSpin,block[i]->linkStrength[j]);
            tot_d_onsiteEnergy+=source_anisotropy;
            if(block[i]->linkedOrb[j]->inBlock>0){
                // target anisotropy
                tot_d_onsiteEnergy+=block[i]->linkedOrb[j]->sDotN*(*p_diagonalDot)(block[i]->perpenSpin,*refDirection,block[i]->linkStrength[j]);
            }else{
                tot_d_onsiteEnergy+=source_anisotropy;
            }
        }
    }
    // single-ion anisotropy
    for(i=0;i<*p_blockLen;i++) tot_d_onsiteEnergy+=getDeltaOnsiteEnergy(block[i]);
    // process the onsite anisotropy
    if(tot_d_onsiteEnergy<=0 || exp(-tot_d_onsiteEnergy)>rand()/(double) RAND_MAX){
        for(i=0;i<*p_blockLen;i++){
            cTimes(&block[i]->transSpin,2);
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
    // clean
    free(refDirection);free(block);free(buffer);
}

void localUpdate(int totOrbs, Orb lattice[], double *p_energy, Vec *p_totSpin){
    //printf("start local updating\n");
    int seedID=(rand()*RAND_MAX+rand())%totOrbs;  // chose one orb Note that WE CANNOT CHOOSE GHOST SPIN, since it's not compatible with local statistics
    //int seedID=totOrbs;
    //printf("considering %d orb, its spin is %.3f %.3f\n",
    //         seedID,lattice[seedID].spin.x,lattice[seedID].spin.y);
    Vec *refDirection=generateRandomVec(); // chose new direction
    //printf("try new spin direction, ref: %.3f %.3f\n",
    //       refDirection->x,refDirection->y);
    double s1n=-2*dot(lattice[seedID].spin,*refDirection);
    //printf("projection s1n: %.3f\n",s1n);
    equal(&lattice[seedID].transSpin,*refDirection);
    cTimes(&lattice[seedID].transSpin,s1n);
    double corr=getDeltaCorrEnergy(lattice+seedID);
    corr+=getDeltaOnsiteEnergy(lattice+seedID);
    
    //printf("lead to the translation spin vector: %.3f %.3f and delta Ecorr: %.3f, transition possibility %.3f P\n",
    //      lattice[seedID].transSpin.x,lattice[seedID].transSpin.y,corr,100*exp(-corr));
    
    if(corr<=0 || exp(-corr)>rand()/(double) RAND_MAX){  // new direction is energertically favoured thus accept directly
        plusEqual(p_totSpin,lattice[seedID].transSpin);
        plusEqual(&lattice[seedID].spin,lattice[seedID].transSpin);
        *p_energy+=corr;
        //printf("Luckily we accept the transition, energy: %.3f\n",*p_energy);
    }
    free(refDirection);
    return;
}

void (*p_mcUpdate)(int totOrbs, Orb lattice[], double*p_energy, Vec *p_totSpin);

PyObject * MCMainFunction(PyObject* self, PyObject* args){
    // read in all parameters
    PyObject* py_algorithm;
    PyObject* py_initSpin;
    PyObject* py_initD;
    PyObject* py_nthermal;
    PyObject* py_nsweep;
    PyObject* py_maxNLinking;
    PyObject* py_ninterval;
    PyObject* py_nlink;
    PyObject* py_linkStrength;
    PyObject* py_linkedOrb;
    PyObject* py_localCircuits;
    PyObject* py_corrOrbPair;
    PyObject* py_nOrbGroup;
    PyObject* py_maxOrbGroupSize;
    PyObject* py_orbGroupList;
    PyObject* py_flunc; 
    PyObject* py_h;
    PyObject* py_rOrb;
    PyObject* py_rOrbCluster;
    PyObject* py_linkedOrb_rnorm;
    PyObject* py_spinFrame;
    PyObject* py_ignoreNonDiagonalJ;
    PyObject* callback;  // callback function
    printf("start parsing...\n");
    PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOOOOOOOOO",
                    &py_algorithm,&py_initSpin,&py_initD,&py_nthermal,&py_nsweep,&py_ninterval,
                    &py_maxNLinking,&py_nlink,&py_linkStrength,&py_linkedOrb,&py_localCircuits,
                    &py_corrOrbPair,
                    &py_nOrbGroup,&py_maxOrbGroupSize,&py_orbGroupList,
                    &py_flunc,&py_h,
                    &py_rOrb,&py_rOrbCluster,&py_linkedOrb_rnorm,
                    &py_spinFrame,&py_ignoreNonDiagonalJ,
                    &callback);

    int algorithm=(int)PyLong_AsLong(py_algorithm); 
    //printf("%d\n",algorithm);
    int nthermal=(int)PyLong_AsLong(py_nthermal);
    int nsweep=(int)PyLong_AsLong(py_nsweep);
    int maxNLinking=(int)PyLong_AsLong(py_maxNLinking);
    long ninterval=PyLong_AsLong(py_ninterval);
    int spinFrame=(int)PyLong_AsLong(py_spinFrame);
    int ignoreNonDiagonalJ=(int)PyLong_AsLong(py_ignoreNonDiagonalJ);

    int totOrbs=(int)PyTuple_Size(py_initSpin);
    double *initSpin=(double*)malloc(totOrbs*sizeof(double));
    for(int iorb=0;iorb<totOrbs;iorb++)initSpin[iorb]=PyFloat_AsDouble(PyTuple_GetItem(py_initSpin,iorb));

    double *initD=(double*)malloc(totOrbs*3*sizeof(double));
    for(int iorb=0;iorb<totOrbs;iorb++)for(int idir=0;idir<3;idir++)
        initD[iorb*3+idir]=PyFloat_AsDouble(PyTuple_GetItem(py_initD,iorb*3+idir));
    
    
    int *nlink=(int*)malloc(totOrbs*sizeof(int));
    double *linkStrength=(double*)malloc(totOrbs*maxNLinking*9*sizeof(double));
    int *linkedOrb=(int*)malloc(totOrbs*maxNLinking*sizeof(int));
    for(int iorb=0;iorb<totOrbs;iorb++){
        nlink[iorb]=(int)PyLong_AsLong(PyTuple_GetItem(py_nlink,iorb));
        for(int ilink=0;ilink<maxNLinking;ilink++){
            linkedOrb[iorb*maxNLinking+ilink]=(int)PyLong_AsLong(PyTuple_GetItem(py_linkedOrb,iorb*maxNLinking+ilink));
            for(int icomp=0;icomp<9;icomp++)
            linkStrength[iorb*maxNLinking*9+ilink*9+icomp]=PyFloat_AsDouble(PyTuple_GetItem(py_linkStrength,iorb*maxNLinking*9+ilink*9+icomp));
        }
    }
    //printf("totOrbs=%d s0=%f D0x=%f nther=%d nst=%d tau=%d maxLink=%d link0=%d J0x=%f\n",totOrbs,initSpin[0],initD[0],nthermal,nsweep,ninterval,maxNLinking,nlink[0],linkStrength[0][0][0]);
    //printf("spinFrame: %d only diagonal: %d\n",spinFrame,ignoreNonDiagonalJ);
    //for(int iorb=0;iorb<nlink[0];iorb++)printf("orb0-orb%d\n",linkedOrb[0][iorb]);
    
    unsigned long long nLocalCircuits=(unsigned long long)PyTuple_Size(py_localCircuits)/3;
    unsigned long long minimalLocalCircuits=1;if(nLocalCircuits>minimalLocalCircuits)minimalLocalCircuits=nLocalCircuits;
    int *localCircuits=(int*)malloc(minimalLocalCircuits*3*sizeof(int));
    for(unsigned long long icircuit=0;icircuit<nLocalCircuits;icircuit++)for(int icomp=0;icomp<3;icomp++)
        localCircuits[icircuit*3+icomp]=(int)PyLong_AsLong(PyTuple_GetItem(py_localCircuits,icircuit*3+icomp));
    //printf("num. of local circuit for topo: %d\n",nLocalCircuits);

    unsigned long long nLat=(unsigned long long)PyTuple_Size(py_corrOrbPair)/2;
    int *corrOrbPair=(int*)malloc(nLat*2*sizeof(int));
    for(unsigned long long ilat=0;ilat<nLat;ilat++){
        for(int icomp=0;icomp<2;icomp++)
        corrOrbPair[ilat*2+icomp]=(int)PyLong_AsLong(PyTuple_GetItem(py_corrOrbPair,ilat*2+icomp));
        //printf("pair%d orb%d-orb%d\n",ilat,corrOrbPair[ilat*2+0],corrOrbPair[ilat*2+1]);
    }
    
    unsigned long long nOrbGroup=(unsigned long long)PyLong_AsLong(py_nOrbGroup);
    unsigned long long maxOrbGroupSize=(unsigned long long)PyLong_AsLong(py_maxOrbGroupSize);
    //printf("nOrbGroup=%d maxOrbGroupSize=%d\n",nOrbGroup,maxOrbGroupSize);
    int *orbGroupList=(int*)malloc(nOrbGroup*maxOrbGroupSize*sizeof(int));
    for(unsigned long long iorbGroup=0;iorbGroup<nOrbGroup;iorbGroup++)for(unsigned long long iorb=0;iorb<maxOrbGroupSize;iorb++)
        orbGroupList[iorbGroup*maxOrbGroupSize+iorb]=(int)PyLong_AsLong(PyTuple_GetItem(py_orbGroupList,iorbGroup*maxOrbGroupSize+iorb));
    
    double flunc = PyFloat_AsDouble(py_flunc);
    double h = PyFloat_AsDouble(py_h);
    //printf("flunc=%f h=%f\n",flunc,h);

    int totOrb_rnorm=(int)PyTuple_Size(py_rOrb);
    int nOrbInCluster=(int)PyTuple_Size(py_rOrbCluster)/totOrb_rnorm;
    //printf("totOrb renorm=%d grain size=%d\n",totOrb_rnorm,nOrbInCluster);
    int *rOrb=(int*)malloc(totOrb_rnorm*sizeof(int));
    int *rOrbCluster=(int*)malloc(totOrb_rnorm*nOrbInCluster*sizeof(int));
    int *linkedOrb_rnorm=(int*)malloc(totOrb_rnorm*maxNLinking*sizeof(int));
    for(int iorb=0;iorb<totOrb_rnorm;iorb++){
        rOrb[iorb]=(int)PyLong_AsLong(PyTuple_GetItem(py_rOrb,iorb));
        for(int iorb_cluster=0;iorb_cluster<nOrbInCluster;iorb_cluster++)
            rOrbCluster[iorb*nOrbInCluster+iorb_cluster]=(int)PyLong_AsLong(PyTuple_GetItem(py_rOrbCluster,iorb*nOrbInCluster+iorb_cluster));
        for(int ilink=0;ilink<maxNLinking;ilink++)
            linkedOrb_rnorm[iorb*maxNLinking+ilink]=(int)PyLong_AsLong(PyTuple_GetItem(py_linkedOrb_rnorm,iorb*maxNLinking+ilink));
    }

    printf("Args parsing success!\n");
    // set algorithm
    p_mcUpdate=localUpdate;
    if(algorithm==1) p_mcUpdate=blockUpdate;

    // set diagonal dot function
    p_diagonalDot=diagonalDot;
    if (ignoreNonDiagonalJ>0) p_diagonalDot=diagonalDot_simple;

    // initialize lattice add one ghost spin for mimicing external field
    Orb *lattice=(Orb*)malloc(totOrbs*sizeof(Orb));
    //printf("hello here is C lib\n");
    establishLattice(lattice, totOrbs, initSpin, initD, flunc, maxNLinking, nlink, linkStrength, h, totOrb_rnorm, nOrbInCluster, rOrb, rOrbCluster);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb, totOrb_rnorm, rOrb, linkedOrb_rnorm);

    // initialize measurement
    double energy=0;
    double *p_energy=&energy;
    for(int i=0;i<totOrbs;i++)*p_energy+=getCorrEnergy(lattice+i);
    *p_energy/=2; // double counting
    for(int i=0;i<totOrbs;i++)*p_energy+=getOnsiteEnergy(lattice+i);

    Vec totSpin;
    totSpin.x=0;totSpin.y=0;
    Vec*p_totSpin=&totSpin;
    for(int i=0;i<totOrbs;i++) plusEqual(p_totSpin, lattice[i].spin);
    
    for(unsigned long long i=0;i<nthermal*ninterval;i++) (*p_mcUpdate)(totOrbs, lattice, p_energy, p_totSpin); //thermalization
    
    printf("start sweeping\n");
    Vec spin_i, spin_i_r;
    Vec spin_j, spin_j_r;
    spin_i.x=0;spin_i.y=0;spin_i_r.x=0;spin_i_r.y=0;
    spin_j.x=0;spin_j.y=0;spin_j_r.x=0;spin_j_r.y=0;
    double spin_ij, spin_ij_r;
    spin_ij=0;spin_ij_r=0;

    double totEnergy=0, totEnergy_r=0;
    double E2=0, E2_r=0;
    double M=0,M2=0,M4=0;
    double M_tmp=0,MdotM_tmp=0,M_tot=0;

    double topological_q=0; // total vortice number

    Vec spin_direction;
    double spin_i_z, spin_j_z, spin_tot_z;
    double spin_i_h, spin_j_h, spin_tot_h;
    spin_i_z=0;spin_j_z=0;spin_tot_z=0;
    spin_i_h=0;spin_j_h=0;spin_tot_h=0;

    // prepare for output spin frame
    int output_per_sweep=nsweep;
    int iFrame=0;
    PyObject *spinFrameData;
    if(spinFrame>0){
        output_per_sweep=nsweep/spinFrame;
        spinFrameData=PyTuple_New(spinFrame);
    }else{
        spinFrameData=PyFloat_FromDouble(0.0);
    }

    // prepare for orb group statistics
    double spinDotSpinBetweenGroup[(nOrbGroup+1)*(nOrbGroup+1)];
    for(unsigned long long i=0;i<(nOrbGroup+1)*(nOrbGroup+1);i++) spinDotSpinBetweenGroup[i]=0.0;
    double spin4Order[nOrbGroup+1];
    for(unsigned long long i=0;i<nOrbGroup+1;i++) spin4Order[i]=0.0;

    for(int i=0;i<nsweep;i++){
        for(unsigned long long j=0;j<ninterval;j++) (*p_mcUpdate)(totOrbs, lattice, p_energy, p_totSpin);
        // record the spin vector field distribution
        if((spinFrame>0) & (i%output_per_sweep==0)){
            PyObject *spinDistribution=PyTuple_New(totOrbs);
            for(int j=0;j<totOrbs;j++){
                PyObject *spinJVec=PyTuple_New(3);
                PyTuple_SetItem(spinJVec, 0, PyFloat_FromDouble(lattice[j].spin.x));
                PyTuple_SetItem(spinJVec, 1, PyFloat_FromDouble(lattice[j].spin.y));
                PyTuple_SetItem(spinJVec, 2, PyFloat_FromDouble(0));
                PyTuple_SetItem(spinDistribution, j, spinJVec);
            }
            PyTuple_SetItem(spinFrameData, iFrame, spinDistribution);
            iFrame+=1;
        }

        // find the main axis
        spin_direction.x=p_totSpin->x;
        spin_direction.y=p_totSpin->y;
        normalize(&spin_direction);

        // spin statistics over space in each frame
        Vec spin_i_avg;
        Vec spin_j_avg;
        spin_i_avg.x=0;
        spin_i_avg.y=0;
        spin_j_avg.x=0;
        spin_j_avg.y=0;
        double spin_ij_avg=0.0;

        double spin_i_z_avg, spin_j_z_avg;
        double spin_i_h_avg, spin_j_h_avg;
        spin_i_z_avg=0;spin_j_z_avg=0;
        spin_i_h_avg=0;spin_j_h_avg=0;
        for(unsigned long long j=0;j<nLat;j++){
            plusEqual(&spin_i_avg, lattice[corrOrbPair[j*2+0]].spin);
            plusEqual(&spin_j_avg, lattice[corrOrbPair[j*2+1]].spin);
            spin_ij_avg+=dot(lattice[corrOrbPair[j*2+0]].spin,lattice[corrOrbPair[j*2+1]].spin);

            // spin along main axis
            spin_i_z_avg+=dot(spin_direction,lattice[corrOrbPair[j*2+0]].spin);
            spin_j_z_avg+=dot(spin_direction,lattice[corrOrbPair[j*2+1]].spin);

            // spin projected to z axis
            spin_i_h_avg+=lattice[corrOrbPair[j*2+0]].spin.x;
            spin_j_h_avg+=lattice[corrOrbPair[j*2+1]].spin.x;
            
        }

        double topological_q_local=0; // xy-model posess no topological charge
        topological_q+=topological_q_local;

        spin_i_z+=spin_i_z_avg/nLat;
        spin_j_z+=spin_j_z_avg/nLat;
        spin_tot_z+=(dot(spin_direction,*p_totSpin)/nLat);
        //if(h<0.00001){// avoid faults time reversal symmetry
        //    spin_i_h+=fabs(spin_i_h_avg)/nLat;
        //    spin_j_h+=fabs(spin_j_h_avg)/nLat;
        //    spin_tot_h+=fabs(p_totSpin->x/nLat);
        //}else{
            spin_i_h+=spin_i_h_avg/nLat;
            spin_j_h+=spin_j_h_avg/nLat;
            spin_tot_h+=p_totSpin->x/nLat;
        //}

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

        double e_avg=*p_energy/totOrbs;
        //printf("e_avg=%.5f\n",e_avg);
        totEnergy+=e_avg;
        E2+=e_avg*e_avg;

        // *************** statistics on renormalized lattice
        // spin statistics over space in each frame
        Vec spin_i_r_avg, spin_j_r_avg;
        spin_i_r_avg.x=0;spin_i_r_avg.y=0;
        spin_j_r_avg.x=0;spin_j_r_avg.y=0;
        int chosed_spin_i=0;
        int chosed_spin_j=0;
        int chosed_spin_ij=0;
        double spin_ij_r_avg=0.0;
        for(unsigned long long j=0;j<nLat;j++){
            Orb *_orb_i=lattice+corrOrbPair[j*2+0];
            Orb *_orb_j=lattice+corrOrbPair[j*2+1];
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

        // orb group statistics
        Vec summedSpinOfGroup[nOrbGroup+1];
        unsigned long long iOrbGroup, jOrbGroup;
        for(iOrbGroup=0;iOrbGroup<nOrbGroup;iOrbGroup++){
            Vec summedSpin;
            summedSpin.x=0;summedSpin.y=0;
            for(unsigned long long k=0;k<maxOrbGroupSize;k++){
                if(orbGroupList[iOrbGroup*maxOrbGroupSize+k]<0) break;
                plusEqual(&summedSpin,lattice[orbGroupList[iOrbGroup*maxOrbGroupSize+k]].spin);
            }
            summedSpinOfGroup[iOrbGroup].x=summedSpin.x;
            summedSpinOfGroup[iOrbGroup].y=summedSpin.y;
        }
        summedSpinOfGroup[iOrbGroup].x=p_totSpin->x/nLat;
        summedSpinOfGroup[iOrbGroup].y=p_totSpin->y/nLat;
        //printf("nOrbGroup: %d\n",nOrbGroup);
        for(iOrbGroup=0;iOrbGroup<(nOrbGroup+1);iOrbGroup++){
            for(jOrbGroup=0;jOrbGroup<(nOrbGroup+1);jOrbGroup++){
                double spin_i_dot_spin_j=dot(summedSpinOfGroup[iOrbGroup],summedSpinOfGroup[jOrbGroup]);
                spinDotSpinBetweenGroup[iOrbGroup*(nOrbGroup+1)+jOrbGroup]+=spin_i_dot_spin_j;
                //printf("i=%d, j=%d, spin_i_dot_j=%.6f\n",iOrbGroup,jOrbGroup,spin_i_dot_spin_j);
                if(iOrbGroup==jOrbGroup) spin4Order[iOrbGroup]+=spin_i_dot_spin_j*spin_i_dot_spin_j;
            }
        }
    }
    //printf("%.3f %.3f\n",spin_i_r.x,spin_i_r.y);
    double U4=(M2/nsweep)*(M2/nsweep)/(M4/nsweep);
    double autoCorr=(MdotM_tmp/nsweep-(M_tot/nsweep)*(M_tot/nsweep));
    PyObject *spinDotSpinBetweenGroup_Tuple;
    spinDotSpinBetweenGroup_Tuple=PyTuple_New((nOrbGroup+2)*(nOrbGroup+1));
    for(int iOrbGroup=0;iOrbGroup<(nOrbGroup+1);iOrbGroup++){
        for(int jOrbGroup=0;jOrbGroup<(nOrbGroup+1);jOrbGroup++){
            int index=iOrbGroup*(nOrbGroup+1)+jOrbGroup;
            double result=spinDotSpinBetweenGroup[index]/nsweep;
            PyTuple_SetItem(spinDotSpinBetweenGroup_Tuple, index, PyFloat_FromDouble(result));
        }
    }
    for(int iOrbGroup=0;iOrbGroup<(nOrbGroup+1);iOrbGroup++){
        int index=(nOrbGroup+1)*(nOrbGroup+1)+iOrbGroup;
        double result=spin4Order[iOrbGroup]/nsweep;
        PyTuple_SetItem(spinDotSpinBetweenGroup_Tuple, index, PyFloat_FromDouble(result));
    }
    if(nOrbGroup==0) spinDotSpinBetweenGroup_Tuple=PyFloat_FromDouble(0);
    PyObject *Data;
    Data=PyTuple_New(29);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i.x/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_i.y/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(0));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(spin_j.x/nsweep));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(spin_j.y/nsweep));
    PyTuple_SetItem(Data, 5, PyFloat_FromDouble(0));
    PyTuple_SetItem(Data, 6, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 7, PyFloat_FromDouble(autoCorr));
    PyTuple_SetItem(Data, 8, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 9, PyFloat_FromDouble(E2/nsweep));
    PyTuple_SetItem(Data,10, PyFloat_FromDouble(U4));
    PyTuple_SetItem(Data,11, PyFloat_FromDouble(spin_i_r.x/nsweep));
    PyTuple_SetItem(Data,12, PyFloat_FromDouble(spin_i_r.y/nsweep));
    PyTuple_SetItem(Data,13, PyFloat_FromDouble(0));
    PyTuple_SetItem(Data,14, PyFloat_FromDouble(spin_j_r.x/nsweep));
    PyTuple_SetItem(Data,15, PyFloat_FromDouble(spin_j_r.y/nsweep));
    PyTuple_SetItem(Data,16, PyFloat_FromDouble(0));
    PyTuple_SetItem(Data,17, PyFloat_FromDouble(spin_ij_r/nsweep));
    PyTuple_SetItem(Data,18, PyFloat_FromDouble(totEnergy_r/nsweep));
    PyTuple_SetItem(Data,19, PyFloat_FromDouble(E2_r/nsweep));
    PyTuple_SetItem(Data,20, PyFloat_FromDouble(spin_i_z/nsweep));
    PyTuple_SetItem(Data,21, PyFloat_FromDouble(spin_j_z/nsweep));
    PyTuple_SetItem(Data,22, PyFloat_FromDouble(spin_tot_z/nsweep));
    PyTuple_SetItem(Data,23, PyFloat_FromDouble(spin_i_h/nsweep));
    PyTuple_SetItem(Data,24, PyFloat_FromDouble(spin_j_h/nsweep));
    PyTuple_SetItem(Data,25, PyFloat_FromDouble(spin_tot_h/nsweep));
    PyTuple_SetItem(Data,26, PyFloat_FromDouble(topological_q/nsweep));
    PyTuple_SetItem(Data, 27, spinFrameData);
    PyTuple_SetItem(Data, 28, spinDotSpinBetweenGroup_Tuple);
    return Data;
}

static PyMethodDef module_methods[] = {
    {"MCMainFunction", MCMainFunction, METH_VARARGS, "the only function in our c lib"},
    {NULL, NULL, 0, NULL} // neccessary to tell python compiler stop here
};

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "xylib",
    "Used to execute the Monte Carlo simulations of xy-model, that is, the O(2) spin model.",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_xylib(void) {
    printf("Initializing xylib...\n");
    PyObject* m;
    m= PyModule_Create(&moduledef);
    return m;
}