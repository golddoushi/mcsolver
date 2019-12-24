#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Orb
{
    int id;
    float spin;
    int nlink;
    float *linkStrength;
    int inBlock;
    struct Orb **linkedOrb;
}Orb;

void establishLattice(Orb *lattice, int totOrbs, float initSpin[totOrbs], int maxNLinking, int nlink[maxNLinking], float linkStrength[totOrbs][maxNLinking]){
    for(int i=0;i<totOrbs;i++){
        lattice[i].id=i;
        lattice[i].spin=initSpin[i];
        lattice[i].nlink=nlink[i];
        lattice[i].linkStrength=linkStrength[i];
    }
}

void establishLinking(Orb *lattice, int totOrbs, int maxNLinking, int nlink[maxNLinking], int linkedOrb[totOrbs][maxNLinking]){
    for(int iorb=0;iorb<totOrbs;iorb++){
        lattice[iorb].linkedOrb=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb][ilink];
        }
    }
}

float getCorrEnergy(Orb *source){
    float corr=0;
    for(int i=0;i<source->nlink;i++){
        corr+=(source->linkStrength[i])*(source->spin)*(source->linkedOrb[i]->spin);
    }
    return corr;
}

int expandBlock(int*beginIndex, int*endIndex, Orb *buffer[], int*blockLen, Orb *block[]){
    //printf("  Buffer: now start and end pt is %d, %d.\n",*beginIndex, *endIndex);
    if(*beginIndex>*endIndex) return 0;

    // FIFO
    Orb *outOrb=buffer[*beginIndex];
    *beginIndex+=1; // pop out the first element
    
    //FILO
    //Orb *outOrb=buffer[*endIndex];
    //*endIndex-=1; // pop out the last element

    int i;
    for(i=0;i<outOrb->nlink;i++){
        Orb *linkedOrb=outOrb->linkedOrb[i];
        //printf("      considering the %d orb which is linking to %d orb, it is %d in block \n", linkedOrb->id, outOrb->id, linkedOrb->inBlock);
        if(linkedOrb->inBlock==0){
            float corr=(outOrb->linkStrength[i])*(outOrb->spin)*(linkedOrb->spin); // bond strength
            //printf("          since it is not in block thus we calc. the correlation energy is %f\n",corr);
            if(corr<0 && (1-exp(2*corr))>rand()/32767.0){
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

void blockUpdate(int totOrbs, Orb lattice[], float*p_energy, float*p_totSpin){
    //printf("one block update step is initializaing...\n");
    Orb *block[totOrbs];
    Orb *buffer[totOrbs];
    int seedID=rand()%totOrbs;
    block[0]=lattice+seedID;
    buffer[0]=lattice+seedID;
    block[0]->inBlock=1;
    int beginIndex=0, endIndex=0, blockLen=1, i;
    int *p_beginIndex=&beginIndex, *p_endIndex=&endIndex, *p_blockLen=&blockLen;

    //printf("the seed Orb is %d\n",block[0]->id);
    while (expandBlock(p_beginIndex, p_endIndex, buffer, p_blockLen, block)==1)
    {
        //printf("    Block size is %d\n",*p_blockLen);
    }
    for(i=0;i<*p_blockLen;i++){
        block[i]->spin*=-1;
        block[i]->inBlock=0;
        *p_totSpin+=(block[i]->spin*2);
    }
    *p_energy=0.;
    for(i=0;i<totOrbs;i++){
        *p_energy+=getCorrEnergy(lattice+i);
    }
    *p_energy/=2;
}

void localUpdate(int totOrbs, Orb lattice[], float *p_energy, float *p_totSpin){
    int seedID=rand()%totOrbs;
    float corr=2*getCorrEnergy(lattice+seedID);

    if(corr>=0){
        lattice[seedID].spin*=-1;
        *p_totSpin+=(lattice[seedID].spin*2);
        *p_energy-=corr;
    }else if (exp(corr)>rand()/32767.0){
        lattice[seedID].spin*=-1;
        *p_totSpin+=(lattice[seedID].spin*2);
        *p_energy-=corr;
    }
}

PyObject * blockUpdateMC(int totOrbs, float initSpin[totOrbs], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], float linkStrength[totOrbs][maxNLinking], int linkedOrb[totOrbs][maxNLinking]){
    //printf("hello block algorithm!\n");
    // initialize lattice
    Orb lattice[totOrbs];
    establishLattice(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb);

    // initialize measurement
    float energy=0, totSpin=0;
    float *p_energy=&energy, *p_totSpin=&totSpin;
    for(int i=0;i<totOrbs;i++) totSpin+=initSpin[i];
    
    // initialize block
    for(int i=0;i<totOrbs;i++) lattice[i].inBlock=0;

    for(int i=0;i<nthermal;i++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin); //thermalization

    // printf("start sweeping\n");
    PyObject *spinData, *energyData;
    spinData=PyTuple_New(nsweep);
    energyData=PyTuple_New(nsweep);
    for(int i=0;i<nsweep;i++){
        blockUpdate(totOrbs, lattice, p_energy, p_totSpin);
        PyTuple_SetItem(spinData, i, PyFloat_FromDouble(*p_totSpin));
        PyTuple_SetItem(energyData, i, PyFloat_FromDouble(*p_energy));
    }
    PyObject *Data;
    Data=PyTuple_New(2);
    PyTuple_SetItem(Data, 0, spinData);
    PyTuple_SetItem(Data, 1, energyData);
    return Data;
}

PyObject * localUpdateMC(int totOrbs, float initSpin[totOrbs], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], float linkStrength[totOrbs][maxNLinking], int linkedOrb[totOrbs][maxNLinking]){
    // initialize lattice
    Orb lattice[totOrbs];
    establishLattice(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb);

    // initialize measurement
    float energy=0, totSpin=0;
    float *p_energy=&energy, *p_totSpin=&totSpin;
    for(int i=0;i<totOrbs;i++) totSpin+=initSpin[i];
    
    for(int i=0;i<totOrbs*nthermal;i++) localUpdate(totOrbs, lattice, p_energy, p_totSpin); //thermalization

    PyObject *spinData, *energyData;
    spinData=PyTuple_New(nsweep);
    energyData=PyTuple_New(nsweep);
    for(int i=0;i<nsweep;i++){
        float energyAvg=0.0;
        float spinAvg=0.0;
        for(int j=0;j<totOrbs;j++){
            localUpdate(totOrbs, lattice, p_energy, p_totSpin);
            energyAvg+=*p_energy;
            spinAvg+=*p_totSpin;
        }
        PyTuple_SetItem(spinData, i, PyFloat_FromDouble(spinAvg/totOrbs));
        PyTuple_SetItem(energyData, i, PyFloat_FromDouble(energyAvg/totOrbs));
    }
    
    PyObject *Data;
    Data=PyTuple_New(2);
    PyTuple_SetItem(Data, 0, spinData);
    PyTuple_SetItem(Data, 1, energyData);
    return Data;
}