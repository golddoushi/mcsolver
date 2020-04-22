#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Orb
{
    int id;
    double spin;
    int nlink;
    double *linkStrength;
    int inBlock;
    struct Orb **linkedOrb;
}Orb;

void establishLattice(Orb *lattice, int totOrbs, double initSpin[totOrbs], int maxNLinking, int nlink[maxNLinking], double linkStrength[totOrbs][maxNLinking]){
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

void establishLatticeWithGhost(Orb *lattice, int totOrbs, double initSpin[totOrbs], int maxNLinking, int nlink[maxNLinking], double linkStrength[totOrbs][maxNLinking], double h){
    for(int i=0;i<totOrbs-1;i++){
        lattice[i].id=i;
        lattice[i].spin=initSpin[i];
        lattice[i].nlink=nlink[i]+1;
        lattice[i].linkStrength=(double*)malloc(lattice[i].nlink*sizeof(double));
        for(int j=0;j<lattice[i].nlink-1;j++)lattice[i].linkStrength[j]=linkStrength[i][j];
        lattice[i].linkStrength[lattice[i].nlink-1]=1;
    }
    lattice[totOrbs-1].id=totOrbs-1;
    lattice[totOrbs-1].spin=h;
    lattice[totOrbs-1].nlink=totOrbs-1;
    lattice[totOrbs-1].linkStrength=(double*)malloc((totOrbs-1)*sizeof(double));
    for(int i=0;i<totOrbs-1;i++) lattice[totOrbs-1].linkStrength[i]=1;
}

void establishLinkingWithGhost(Orb *lattice, int totOrbs, int maxNLinking, int nlink[maxNLinking], int linkedOrb[totOrbs][maxNLinking]){
    for(int iorb=0;iorb<totOrbs-1;iorb++){
        lattice[iorb].linkedOrb=(Orb**)malloc(lattice[iorb].nlink*sizeof(Orb*));
        for(int ilink=0;ilink<lattice[iorb].nlink-1;ilink++){
            lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb][ilink];
        }
        lattice[iorb].linkedOrb[lattice[iorb].nlink-1]=lattice+totOrbs-1;
    }
    lattice[totOrbs-1].linkedOrb=(Orb**)malloc((totOrbs-1)*sizeof(Orb*));
    for(int i=0;i<totOrbs-1;i++) lattice[totOrbs-1].linkedOrb[i]=lattice+i;
}

double getCorrEnergy(Orb *source){
    double corr=0;
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
    //printf("there are %d linked orbs\n",outOrb->nlink);
    for(i=0;i<outOrb->nlink;i++){
        Orb *linkedOrb=outOrb->linkedOrb[i];
        //printf("      considering the %d orb which is linking to %d orb, it is %d in block \n", linkedOrb->id, outOrb->id, linkedOrb->inBlock);
        if(linkedOrb->inBlock==0){
            double corr=(outOrb->linkStrength[i])*(outOrb->spin)*(linkedOrb->spin); // bond strength
            //printf("          since it is not in block thus we calc. the correlation energy is %f\n",corr);
            if(corr<0 && (1-exp(2*corr))>rand()/32767.0){
                //printf("          -->>fortunately it is added to block, Padd=%f\n",(1-exp(2*corr)));
                // update block
                *blockLen+=1;
                block[*blockLen-1]=linkedOrb;
                linkedOrb->inBlock=1;  // register into block
                // update buffer
                *endIndex+=1;
                buffer[*endIndex]=linkedOrb;
            }//else{
            //    printf("          -->>unfortunately it is not added to block, Padd=%f\n",(1-exp(2*corr)));
            //}
        }
    }
    return 1;
}

void blockUpdate(int totOrbs, Orb lattice[], double*p_energy, double*p_totSpin){
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
        if(block[i]->id<totOrbs-1) *p_totSpin+=(block[i]->spin*2);
    }
    *p_energy=0.;
    for(i=0;i<totOrbs-1;i++){
        *p_energy+=getCorrEnergy(lattice+i);
    }
    *p_energy/=2;
}

void localUpdate(int totOrbs, Orb lattice[], double *p_energy, double *p_totSpin, double h){
    int seedID=rand()%totOrbs;
    double corr=2*(getCorrEnergy(lattice+seedID)-h*lattice[seedID].spin);

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
 
PyObject * blockUpdateMC(int totOrbs, double initSpin[totOrbs], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking], int linkedOrb[totOrbs][maxNLinking],
                   int ninterval, int nLat, int corrOrbPair[nLat][2], double h){
    //printf("hello block algorithm!\n");
    // initialize lattice including one ghost spin
    totOrbs+=1;
    Orb lattice[totOrbs];
    establishLatticeWithGhost(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength, h);
    establishLinkingWithGhost(lattice, totOrbs, maxNLinking, nlink, linkedOrb);

    // initialize measurement
    double energy=0, totSpin=0;
    double *p_energy=&energy, *p_totSpin=&totSpin;
    for(int i=0;i<totOrbs-1;i++) totSpin+=initSpin[i];
    
    // initialize block
    for(int i=0;i<totOrbs;i++) lattice[i].inBlock=0;
    //printf("initialization success\n");

    for(int i=0;i<nthermal*ninterval;i++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin); //thermalization
    //printf("thermalization finished\n");

    // printf("start sweeping\n");
    //PyObject *spinData, *energyData, *corrData;
    //spinData=PyTuple_New(nsweep);
    //energyData=PyTuple_New(nsweep);
    //corrData=PyTuple_New(nsweep);

    double spin_i=0;
    double spin_j=0;
    double spin_ij=0;
    double totEnergy=0;
    double E2=0;
    for(int i=0;i<nsweep;i++){
        // spin statistics over space in each frame
        double spin_i_avg=0;
        double spin_j_avg=0;
        double corrAvg=0.0;
        for(int j=0;j<nLat;j++){
            double si_tmp=lattice[corrOrbPair[j][0]].spin;
            double sj_tmp=lattice[corrOrbPair[j][1]].spin;
            spin_i_avg+=si_tmp;
            spin_j_avg+=sj_tmp;
            corrAvg+=si_tmp*sj_tmp;
        }

        spin_i+=fabs(spin_i_avg)/nLat;
        spin_j+=fabs(spin_j_avg)/nLat;
        spin_ij+=corrAvg/nLat;

        // energy stored in each frame
        double e_avg=*p_energy/nLat;
        totEnergy+=e_avg;
        E2+=e_avg*e_avg;

        //blockUpdate(totOrbs, lattice, p_energy, p_totSpin);
        for(int j=0;j<ninterval;j++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin);
        //PyTuple_SetItem(spinData, i, PyFloat_FromDouble(*p_totSpin));
        //PyTuple_SetItem(energyData, i, PyFloat_FromDouble(*p_energy));
        //PyTuple_SetItem(corrData, i, PyFloat_FromDouble(corrAvg/nLat));
    }

    PyObject *Data;
    Data=PyTuple_New(5);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_j/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(E2/nsweep));
    return Data;
}

PyObject * localUpdateMC(int totOrbs, double initSpin[totOrbs], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking], int linkedOrb[totOrbs][maxNLinking],
                   int ninterval, int nLat, int corrOrbPair[nLat][2], double h){
    // initialize lattice
    Orb lattice[totOrbs];
    establishLattice(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb);

    // initialize measurement
    double energy=0, totSpin=0;
    double *p_energy=&energy, *p_totSpin=&totSpin;
    for(int i=0;i<ninterval;i++) totSpin+=initSpin[i];
    
    for(int i=0;i<nthermal*ninterval;i++) localUpdate(totOrbs, lattice, p_energy, p_totSpin, h); //thermalization

    //PyObject *spinData, *energyData, *corrData;
    //spinData=PyTuple_New(nsweep);
    //energyData=PyTuple_New(nsweep);
    //corrData=PyTuple_New(nsweep);

    double spin_i=0;
    double spin_j=0;
    double spin_ij=0;
    double totEnergy=0;
    double E2=0;
    for(int i=0;i<nsweep;i++){
        // spin statistics over space in each frame
        double spin_i_avg=0;
        double spin_j_avg=0;
        double corrAvg=0.0;
        for(int j=0;j<nLat;j++){
            double si_tmp=lattice[corrOrbPair[j][0]].spin;
            double sj_tmp=lattice[corrOrbPair[j][1]].spin;
            spin_i_avg+=si_tmp;
            spin_j_avg+=sj_tmp;
            corrAvg+=si_tmp*sj_tmp;
        }
        spin_i+=fabs(spin_i_avg)/nLat;
        spin_j+=fabs(spin_j_avg)/nLat;
        spin_ij+=corrAvg/nLat;

        // energy stored in each frame
        double e_avg=*p_energy/nLat;
        totEnergy+=e_avg;
        E2+=e_avg*e_avg;

        for(int j=0;j<ninterval;j++) localUpdate(totOrbs, lattice, p_energy, p_totSpin, h);
        
        //PyTuple_SetItem(spinData, i, PyFloat_FromDouble(*p_totSpin/ninterval));
        //PyTuple_SetItem(energyData, i, PyFloat_FromDouble(*p_energy/ninterval));
        //PyTuple_SetItem(corrData, i, PyFloat_FromDouble(corrAvg/nLat));
    }
    
    PyObject *Data;
    Data=PyTuple_New(5);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_j/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(E2/nsweep));
    return Data;
}