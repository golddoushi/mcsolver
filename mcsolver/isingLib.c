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

    int chosen;
    struct Orb **linkedOrb_rnorm;
    int nOrbInCluster;
    struct Orb **orb_cluster;
}Orb;

void establishLattice(Orb *lattice, int totOrbs, double initSpin[totOrbs], int maxNLinking, int nlink[maxNLinking], double linkStrength[totOrbs][maxNLinking],
                   int totOrb_rnorm, int nOrbInCluster, int rOrb[totOrb_rnorm], int rOrbCluster[totOrb_rnorm][nOrbInCluster]){
    for(int i=0;i<totOrbs;i++){
        lattice[i].id=i;
        lattice[i].spin=initSpin[i];
        lattice[i].nlink=nlink[i];
        lattice[i].linkStrength=linkStrength[i];
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
}

void establishLinking(Orb *lattice, int totOrbs, int maxNLinking, int nlink[maxNLinking], int linkedOrb[totOrbs][maxNLinking],
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

void establishLatticeWithGhost(Orb *lattice, int totOrbs, double initSpin[totOrbs], int maxNLinking, int nlink[maxNLinking], double linkStrength[totOrbs][maxNLinking], double h,
                   int totOrb_rnorm, int rOrb[totOrb_rnorm]){
    for(int i=0;i<totOrbs-1;i++){
        lattice[i].id=i;
        lattice[i].spin=initSpin[i];
        lattice[i].nlink=nlink[i]+1;
        lattice[i].linkStrength=(double*)malloc(lattice[i].nlink*sizeof(double));
        for(int j=0;j<lattice[i].nlink-1;j++)lattice[i].linkStrength[j]=linkStrength[i][j];
        lattice[i].linkStrength[lattice[i].nlink-1]=1;
        lattice[i].chosen=0;
    }
    for(int i=0;i<totOrb_rnorm;i++){
        lattice[rOrb[i]].chosen=1;
    }
    lattice[totOrbs-1].id=totOrbs-1;
    lattice[totOrbs-1].spin=h;
    lattice[totOrbs-1].nlink=totOrbs-1;
    lattice[totOrbs-1].linkStrength=(double*)malloc((totOrbs-1)*sizeof(double));
    for(int i=0;i<totOrbs-1;i++) lattice[totOrbs-1].linkStrength[i]=1;
}

void establishLinkingWithGhost(Orb *lattice, int totOrbs, int maxNLinking, int nlink[maxNLinking], int linkedOrb[totOrbs][maxNLinking],
                   int totOrb_rnorm, int rOrb[totOrb_rnorm], int linkedOrb_rnorm[totOrb_rnorm][maxNLinking]){
    for(int iorb=0;iorb<totOrbs-1;iorb++){
        lattice[iorb].linkedOrb=(Orb**)malloc(lattice[iorb].nlink*sizeof(Orb*));
        for(int ilink=0;ilink<lattice[iorb].nlink-1;ilink++){
            lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb][ilink];
        }
        lattice[iorb].linkedOrb[lattice[iorb].nlink-1]=lattice+totOrbs-1;
    }
    for(int i=0;i<totOrb_rnorm;i++){
        int iorb=rOrb[i];
        //printf("step %d, In renormalized lattice, orb%d is involved by %d bonds\n",i,lattice[iorb].id, nlink[iorb]);
        lattice[iorb].linkedOrb_rnorm=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb_rnorm[ilink]=lattice+linkedOrb_rnorm[i][ilink];
            //printf("%d\n",linkedOrb_rnorm[i][ilink]);
            //printf("orb%d with strength %.3f\n",lattice[iorb].linkedOrb_rnorm[ilink]->id,lattice[iorb].linkStrength[ilink]);
        }
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

double getMajoritySpin(Orb *_orb){
    //printf("start calc majority spin in block, %d orbs in total, centering orb%d\n",_orb->nOrbInCluster,_orb->id);
    double avg_spin=0.0;
    for(int ispin=0;ispin<_orb->nOrbInCluster;ispin++){
        avg_spin+=_orb->orb_cluster[ispin]->spin;
    }
    if (avg_spin>0){
        return fabs(_orb->spin);
    }else if (avg_spin<0){
        return -fabs(_orb->spin);
    }else{
        if (rand()/(double) RAND_MAX>0.5){
            return fabs(_orb->spin);
        }else{
            return -fabs(_orb->spin);
        }
    }
}

double getCorrEnergy_rnorm(Orb *source){
    //printf("now we are calc. the corr. energy to orb%d\n",source->id);
    double corr=0;
    double avg_spin_source=getMajoritySpin(source);
    for(int i=0;i<source->nlink;i++){
        //printf("link to orb%d\n",source->linkedOrb_rnorm[i]->id);
        double avg_spin_target=getMajoritySpin(source->linkedOrb_rnorm[i]);
        corr+=(source->linkStrength[i])*avg_spin_source*avg_spin_target;
    }
    //printf("Ecorr=%.3f\n",corr);
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
            if(corr<0 && (1-exp(2*corr))>rand()/(double) RAND_MAX){
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

void blockUpdate(int totOrbs, Orb lattice[], double*p_energy, double*p_totSpin, double*p_energy_rnorm){
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
    // energy of renormalized lattice
    //printf("try to evaluate the energy of renormalized system\n");
    *p_energy_rnorm=0.;
    for(i=0;i<totOrbs-1;i++){
        if(lattice[i].chosen>0) *p_energy_rnorm+=getCorrEnergy_rnorm(lattice+i);
    }
    *p_energy_rnorm/=2;
}

void localUpdate(int totOrbs, Orb lattice[], double *p_energy, double *p_totSpin, double h){
    int seedID=rand()%totOrbs;
    double corr=2*(getCorrEnergy(lattice+seedID)-h*lattice[seedID].spin);

    if(corr>=0){
        lattice[seedID].spin*=-1;
        *p_totSpin+=(lattice[seedID].spin*2);
        *p_energy-=corr;
    }else if (exp(corr)>rand()/(double) RAND_MAX){
        lattice[seedID].spin*=-1;
        *p_totSpin+=(lattice[seedID].spin*2);
        *p_energy-=corr;
    }
}
 
PyObject * blockUpdateMC(int totOrbs, double initSpin[totOrbs], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking], int linkedOrb[totOrbs][maxNLinking],
                   int ninterval, int nLat, int corrOrbPair[nLat][2], double h,
                   int totOrb_rnorm, int nOrbInCluster, int rOrb[totOrb_rnorm], int rOrbCluster[totOrb_rnorm][nOrbInCluster], int linkedOrb_rnorm[totOrb_rnorm][maxNLinking]){
    //printf("hello block algorithm!\n");
    // initialize lattice including one ghost spin
    //totOrbs+=1;
    Orb lattice[totOrbs];
    //establishLatticeWithGhost(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength, h, totOrb_rnorm, rOrb);
    //establishLinkingWithGhost(lattice, totOrbs, maxNLinking, nlink, linkedOrb, totOrb_rnorm, rOrb, linkedOrb_rnorm);
    establishLattice(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength, totOrb_rnorm, nOrbInCluster, rOrb, rOrbCluster);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb, totOrb_rnorm, rOrb, linkedOrb_rnorm);

    // initialize measurement
    double energy=0, totSpin=0, energy_rnorm=0;
    double *p_energy=&energy, *p_totSpin=&totSpin, *p_energy_rnorm=&energy_rnorm;
    for(int i=0;i<totOrbs-1;i++) totSpin+=initSpin[i];
    
    // initialize block
    for(int i=0;i<totOrbs;i++) lattice[i].inBlock=0;
    //printf("initialization success\n");

    for(int i=0;i<nthermal*ninterval;i++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin, p_energy_rnorm); //thermalization
    //printf("thermalization finished\n");

    //printf("start sweeping\n");
    double spin_i=0;
    double spin_j=0;
    double spin_ij=0;
    double totEnergy=0, totEnergy_rnorm=0;
    double E2=0, E2_rnorm=0;
    double M=0,M2=0,M4=0;
    double M_tmp=0,MdotM_tmp=0,M_tot=0;

    double spin_tot=0;
    for(int i=0;i<nsweep;i++){
        for(int j=0;j<ninterval;j++) blockUpdate(totOrbs, lattice, p_energy, p_totSpin, p_energy_rnorm); // one sweep
        //printf("sweep%d finished\n",i);
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
        spin_tot+=totSpin;

        M=spin_i_avg/nLat;
        M2+=M*M;
        M4+=M*M*M*M;
        //calc auto-correlation
        M_tot+=M;
        MdotM_tmp+=M_tmp*M;
        M_tmp=M;

        spin_i+=fabs(spin_i_avg)/nLat;
        spin_j+=fabs(spin_j_avg)/nLat;
        spin_ij+=corrAvg/nLat;

        // energy stored in each frame
        //printf("for full lattice with %d orbs, energy is %.3f\n",totOrbs,*p_energy);
        //printf("for renormalized lattice with %d orbs, energy is %.3f\n",totOrb_rnorm,*p_energy_rnorm);
        double e_avg=*p_energy/totOrbs, e_avg_rnorm=*p_energy_rnorm/totOrb_rnorm;
        //printf("avg: %.3f, %.3f\n",e_avg,e_avg_rnorm);
        totEnergy+=e_avg;
        totEnergy_rnorm+=e_avg_rnorm;
        E2+=e_avg*e_avg;
        E2_rnorm+=e_avg_rnorm*e_avg_rnorm;
    }
    //printf("average energy of full lattice: %.3f\n",totEnergy/nsweep);
    //printf("average energy of renormalized lattice: %.3f\n",totEnergy_rnorm/nsweep);
    double U4=(M2/nsweep)*(M2/nsweep)/(M4/nsweep);
    double autoCorr=(MdotM_tmp/nsweep-(M_tot/nsweep)*(M_tot/nsweep));
    PyObject *Data;
    Data=PyTuple_New(10);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_j/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(autoCorr));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 5, PyFloat_FromDouble(E2/nsweep));
    PyTuple_SetItem(Data, 6, PyFloat_FromDouble(totEnergy_rnorm/nsweep));
    PyTuple_SetItem(Data, 7, PyFloat_FromDouble(E2_rnorm/nsweep));
    PyTuple_SetItem(Data, 8, PyFloat_FromDouble(U4));
    PyTuple_SetItem(Data, 9, PyFloat_FromDouble(spin_tot/nsweep/nLat));
    return Data;
}

PyObject * localUpdateMC(int totOrbs, double initSpin[totOrbs], int nthermal, int nsweep, 
                   int maxNLinking, int nlink[totOrbs], double linkStrength[totOrbs][maxNLinking], int linkedOrb[totOrbs][maxNLinking],
                   int ninterval, int nLat, int corrOrbPair[nLat][2], double h,
                   int totOrb_rnorm, int nOrbInCluster, int rOrb[totOrb_rnorm], int rOrbCluster[totOrb_rnorm][nOrbInCluster], int linkedOrb_rnorm[totOrb_rnorm][maxNLinking]){
    // initialize lattice 
    Orb lattice[totOrbs];
    establishLattice(lattice, totOrbs, initSpin, maxNLinking, nlink, linkStrength, totOrb_rnorm, nOrbInCluster, rOrb, rOrbCluster);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb, totOrb_rnorm, rOrb, linkedOrb_rnorm);

    // initialize measurement
    double energy=0, totSpin=0, energy_rnorm=0;;
    double *p_energy=&energy, *p_totSpin=&totSpin, *p_energy_rnorm=&energy_rnorm;
    for(int i=0;i<ninterval;i++) totSpin+=initSpin[i];
    
    for(int i=0;i<nthermal*ninterval;i++) localUpdate(totOrbs, lattice, p_energy, p_totSpin, h); //thermalization

    double spin_i=0;
    double spin_j=0;
    double spin_ij=0;
    double totEnergy=0, totEnergy_rnorm=0;
    double E2=0, E2_rnorm=0;
    *p_energy=0;
    double M=0,M2=0,M4=0;
    double M_tmp=0,MdotM_tmp=0,M_tot=0;

    double spin_tot=0;
    for(int i=0;i<totOrbs;i++){ // calc. energy in renormalized system
        *p_energy+=getCorrEnergy(lattice+i);
    }
    *p_energy/=2;
    for(int i=0;i<nsweep;i++){
        //printf("sweep%d/%d\n",i,nsweep);
        for(int j=0;j<ninterval;j++) localUpdate(totOrbs, lattice, p_energy, p_totSpin, h); // one sweep
        *p_energy_rnorm=0;
        for(int j=0;j<totOrbs;j++){ // calc. energy in renormalized system
            if(lattice[j].chosen>0) *p_energy_rnorm+=getCorrEnergy_rnorm(lattice+j);
        }
        *p_energy_rnorm/=2;
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
        spin_tot+=totSpin;
        M=spin_i_avg/nLat;
        M2+=M*M;
        M4+=M*M*M*M;
        //calc auto-correlation
        M_tot+=M;
        MdotM_tmp+=M_tmp*M;
        M_tmp=M;

        spin_i+=(spin_i_avg)/nLat;//fabs(spin_i_avg)/nLat;
        spin_j+=(spin_j_avg)/nLat;//fabs(spin_j_avg)/nLat;
        spin_ij+=corrAvg/nLat;

        // energy stored in each frame
        //printf("for full lattice with %d orbs, energy is %.3f\n",totOrbs,*p_energy);
        //printf("for renormalized lattice with %d orbs, energy is %.3f\n",totOrb_rnorm,*p_energy_rnorm);
        double e_avg=*p_energy/totOrbs, e_avg_rnorm=*p_energy_rnorm/totOrb_rnorm;
        //printf("avg: %.3f, %.3f\n",e_avg,e_avg_rnorm);
        totEnergy+=e_avg;
        totEnergy_rnorm+=e_avg_rnorm;
        E2+=e_avg*e_avg;
        E2_rnorm+=e_avg_rnorm*e_avg_rnorm;
    }
    double U4=(M2/nsweep)*(M2/nsweep)/(M4/nsweep);
    double autoCorr=(MdotM_tmp/nsweep-(M_tot/nsweep)*(M_tot/nsweep));
    PyObject *Data;
    Data=PyTuple_New(10);
    PyTuple_SetItem(Data, 0, PyFloat_FromDouble(spin_i/nsweep));
    PyTuple_SetItem(Data, 1, PyFloat_FromDouble(spin_j/nsweep));
    PyTuple_SetItem(Data, 2, PyFloat_FromDouble(spin_ij/nsweep));
    PyTuple_SetItem(Data, 3, PyFloat_FromDouble(autoCorr));
    PyTuple_SetItem(Data, 4, PyFloat_FromDouble(totEnergy/nsweep));
    PyTuple_SetItem(Data, 5, PyFloat_FromDouble(E2/nsweep));
    PyTuple_SetItem(Data, 6, PyFloat_FromDouble(totEnergy_rnorm/nsweep));
    PyTuple_SetItem(Data, 7, PyFloat_FromDouble(E2_rnorm/nsweep));
    PyTuple_SetItem(Data, 8, PyFloat_FromDouble(U4));
    PyTuple_SetItem(Data, 9, PyFloat_FromDouble(spin_tot/nsweep/nLat));
    return Data;
}