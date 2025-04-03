#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct Orb
{
    int id;
    double spin;
    double h;
    int nlink;
    double *linkStrength;
    int inBlock;
    struct Orb **linkedOrb;

    int chosen;
    struct Orb **linkedOrb_rnorm;
    int nOrbInCluster;
    struct Orb **orb_cluster;
}Orb;

void establishLattice(Orb *lattice, int totOrbs, double initSpin[totOrbs], double h, int maxNLinking, int nlink[maxNLinking], double linkStrength[totOrbs*maxNLinking],
                   int totOrb_rnorm, int nOrbInCluster, int rOrb[totOrb_rnorm], int rOrbCluster[totOrb_rnorm*nOrbInCluster]){
    for(int i=0;i<totOrbs;i++){
        lattice[i].id=i;
        lattice[i].spin=initSpin[i];
        lattice[i].h=h;
        lattice[i].nlink=nlink[i];
        lattice[i].linkStrength=linkStrength+i;
        lattice[i].chosen=0;
        lattice[i].inBlock=0;
    }
    //printf("verify lattice\n");
    //for(int i=0;i<totOrbs;i++){
    //    printf("orb%d spin: %.3f nlink: %d\n",lattice[i].id,lattice[i].spin,lattice[i].nlink);
    //    for(int j=0;j<nlink[i];j++) printf(" coupling%d: %.6f\n",j,lattice[i].linkStrength[j]);
    //}
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
}

void establishLinking(Orb *lattice, int totOrbs, int maxNLinking, int nlink[maxNLinking], int linkedOrb[totOrbs*maxNLinking],
                   int totOrb_rnorm, int rOrb[totOrb_rnorm], int linkedOrb_rnorm[totOrb_rnorm*maxNLinking]){
    for(int iorb=0;iorb<totOrbs;iorb++){
        lattice[iorb].linkedOrb=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb*maxNLinking+ilink];
        }
    }
    //printf("verify linking\n");
    //for(int iorb=0;iorb<totOrbs;iorb++){
    //    for(int ilink=0;ilink<nlink[iorb];ilink++){
    //        lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb][ilink];
    //        printf("orb%d ~ orb%d strength: %.6f\n",lattice[iorb].id,lattice[iorb].linkedOrb[ilink]->id,lattice[iorb].linkStrength[ilink]);
    //    }
    //}
    for(int i=0;i<totOrb_rnorm;i++){
        int iorb=rOrb[i];
        lattice[iorb].linkedOrb_rnorm=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb_rnorm[ilink]=lattice+linkedOrb_rnorm[i*maxNLinking+ilink];
        }
    }
}

void establishLatticeWithGhost(Orb *lattice, int totOrbs, double initSpin[totOrbs], int maxNLinking, int nlink[maxNLinking], double linkStrength[totOrbs*maxNLinking], double h,
                   int totOrb_rnorm, int rOrb[totOrb_rnorm]){
    for(int i=0;i<totOrbs-1;i++){
        lattice[i].id=i;
        lattice[i].spin=initSpin[i];
        lattice[i].nlink=nlink[i]+1;
        lattice[i].linkStrength=(double*)malloc(lattice[i].nlink*sizeof(double));
        for(int j=0;j<lattice[i].nlink-1;j++)lattice[i].linkStrength[j]=linkStrength[i*maxNLinking+j];
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

void establishLinkingWithGhost(Orb *lattice, int totOrbs, int maxNLinking, int nlink[maxNLinking], int linkedOrb[totOrbs*maxNLinking],
                   int totOrb_rnorm, int rOrb[totOrb_rnorm], int linkedOrb_rnorm[totOrb_rnorm*maxNLinking]){
    for(int iorb=0;iorb<totOrbs-1;iorb++){
        lattice[iorb].linkedOrb=(Orb**)malloc(lattice[iorb].nlink*sizeof(Orb*));
        for(int ilink=0;ilink<lattice[iorb].nlink-1;ilink++){
            lattice[iorb].linkedOrb[ilink]=lattice+linkedOrb[iorb*maxNLinking+ilink];
        }
        lattice[iorb].linkedOrb[lattice[iorb].nlink-1]=lattice+totOrbs-1;
    }
    for(int i=0;i<totOrb_rnorm;i++){
        int iorb=rOrb[i];
        //printf("step %d, In renormalized lattice, orb%d is involved by %d bonds\n",i,lattice[iorb].id, nlink[iorb]);
        lattice[iorb].linkedOrb_rnorm=(Orb**)malloc(nlink[iorb]*sizeof(Orb*));
        for(int ilink=0;ilink<nlink[iorb];ilink++){
            lattice[iorb].linkedOrb_rnorm[ilink]=lattice+linkedOrb_rnorm[i*maxNLinking+ilink];
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

double getDeltaOnsiteEnergy(Orb *source){ // Only field term contributes, 
    return 2*source->h*source->spin;
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

void blockUpdate(int totOrbs, Orb lattice[], double*p_energy, double*p_totSpin){
    //printf("one block update step is initializaing...\n");
    for(int i=0;i<totOrbs;i++) lattice[i].inBlock=0; // initialize all orb status
    Orb **block=(Orb**)malloc(totOrbs*sizeof(Orb*));
    Orb **buffer=(Orb**)malloc(totOrbs*sizeof(Orb*));
    unsigned long long r1=(unsigned long long)rand();
    unsigned long long r2=(unsigned long long)rand();
    unsigned long long seedID=(r1*RAND_MAX+r2)%totOrbs;  // chose one orb
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

    double tot_d_onsiteEnergy=0; // case field on
    // single-ion anisotropy
    for(i=0;i<*p_blockLen;i++) tot_d_onsiteEnergy+=getDeltaOnsiteEnergy(block[i]);
    if(tot_d_onsiteEnergy<=0 || exp(-tot_d_onsiteEnergy)>rand()/(double) RAND_MAX){
        for(i=0;i<*p_blockLen;i++){
            block[i]->spin*=-1;
            if(block[i]->id<totOrbs-1) *p_totSpin+=(block[i]->spin*2);
        }
        *p_energy=0.;
        for(i=0;i<totOrbs;i++){
            *p_energy+=getCorrEnergy(lattice+i)/2-(lattice+i)->h*(lattice+i)->spin;
        }
    }
    free(block);free(buffer);
}

void localUpdate(int totOrbs, Orb lattice[], double *p_energy, double *p_totSpin){
    unsigned long long r1=(unsigned long long)rand();
    unsigned long long r2=(unsigned long long)rand();
    unsigned long long seedID=(r1*RAND_MAX+r2)%totOrbs;  // chose one orb
    double corr=2*(getCorrEnergy(lattice+seedID)-lattice[seedID].h*lattice[seedID].spin);
    //printf("local update: try to flip orb%d, corr=%.6f\n",lattice[seedID].id,corr);
    if(corr>=0){
        lattice[seedID].spin*=-1;
        *p_totSpin+=(lattice[seedID].spin*2);
        *p_energy-=corr;
    }else if (exp(corr)>rand()/(double) RAND_MAX){
        lattice[seedID].spin*=-1;
        *p_totSpin+=(lattice[seedID].spin*2);
        *p_energy-=corr;
    }
    //printf("local update finished\n");
}
 
// interface to block update and local update algorithm
void (*p_mcUpdate)(int totOrbs, Orb lattice[], double *p_energy, double *p_totSpin);

PyObject * MCMainFunction(PyObject* self, PyObject* args){
    // read in all parameters
    PyObject* py_algorithm;
    PyObject* py_initSpin;
    PyObject* py_nthermal;
    PyObject* py_nsweep;
    PyObject* py_maxNLinking;
    PyObject* py_ninterval;
    PyObject* py_nlink;
    PyObject* py_linkStrength;
    PyObject* py_linkedOrb;
    PyObject* py_corrOrbPair;
    PyObject* py_h;
    PyObject* py_rOrb;
    PyObject* py_rOrbCluster;
    PyObject* py_linkedOrb_rnorm;
    PyObject* py_spinFrame;
    PyObject* callback;  // callback function
    PyArg_ParseTuple(args,"OOOOOOOOOOOOOOOO",
                    &py_algorithm,&py_initSpin,&py_nthermal,&py_nsweep,&py_ninterval,
                    &py_maxNLinking,&py_nlink,&py_linkStrength,&py_linkedOrb,
                    &py_corrOrbPair,
                    &py_h,
                    &py_rOrb,&py_rOrbCluster,&py_linkedOrb_rnorm,
                    &py_spinFrame,
                    &callback);
    
    int algorithm=(int)PyLong_AsLong(py_algorithm); 
    //printf("%d\n",algorithm);
    int nthermal=(int)PyLong_AsLong(py_nthermal);
    int nsweep=(int)PyLong_AsLong(py_nsweep);
    int maxNLinking=(int)PyLong_AsLong(py_maxNLinking);
    int ninterval=(int)PyLong_AsLong(py_ninterval);
    int spinFrame=(int)PyLong_AsLong(py_spinFrame);

    int totOrbs=(int)PyTuple_Size(py_initSpin);
    double *initSpin=(double*)malloc(totOrbs*sizeof(double));
    for(int iorb=0;iorb<totOrbs;iorb++)initSpin[iorb]=PyFloat_AsDouble(PyTuple_GetItem(py_initSpin,iorb));
    
    int *nlink=(int*)malloc(totOrbs*sizeof(int));
    double *linkStrength=(double*)malloc(totOrbs*maxNLinking*sizeof(double));
    int *linkedOrb=(int*)malloc(totOrbs*maxNLinking*sizeof(int));
    for(int iorb=0;iorb<totOrbs;iorb++){
        nlink[iorb]=(int)PyLong_AsLong(PyTuple_GetItem(py_nlink,iorb));
        for(int ilink=0;ilink<maxNLinking;ilink++){
            linkedOrb[iorb*maxNLinking+ilink]=(int)PyLong_AsLong(PyTuple_GetItem(py_linkedOrb,iorb*maxNLinking+ilink));
            linkStrength[iorb*maxNLinking+ilink]=PyFloat_AsDouble(PyTuple_GetItem(py_linkStrength,iorb*maxNLinking+ilink));
        }
    }
    //printf("totOrbs=%d s0=%f nther=%d nst=%d tau=%d maxLink=%d link0=%d J=%f\n",totOrbs,initSpin[0],nthermal,nsweep,ninterval,maxNLinking,nlink[0],linkStrength[0][0]);
    //printf("spinFrame: %d\n",spinFrame);
    //for(int iorb=0;iorb<nlink[0];iorb++)printf("orb0-orb%d\n",linkedOrb[0][iorb]);

    int nLat=(int)PyTuple_Size(py_corrOrbPair)/2;
    int *corrOrbPair=(int*)malloc(nLat*2*sizeof(int));
    for(int ilat=0;ilat<nLat;ilat++){
        for(int icomp=0;icomp<2;icomp++)
        corrOrbPair[ilat*2+icomp]=(int)PyLong_AsLong(PyTuple_GetItem(py_corrOrbPair,ilat*2+icomp));
        //printf("pair%d orb%d-orb%d\n",ilat,corrOrbPair[ilat][0],corrOrbPair[ilat][1]);
    }
    
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
    if (algorithm==1) p_mcUpdate=blockUpdate;

    // initialize lattice 
    Orb *lattice=(Orb*)malloc(totOrbs*sizeof(Orb));
    establishLattice(lattice, totOrbs, initSpin, h, maxNLinking, nlink, linkStrength, totOrb_rnorm, nOrbInCluster, rOrb, rOrbCluster);
    establishLinking(lattice, totOrbs, maxNLinking, nlink, linkedOrb, totOrb_rnorm, rOrb, linkedOrb_rnorm);

    // initialize measurement
    double energy=0, totSpin=0, energy_rnorm=0;;
    double *p_energy=&energy, *p_totSpin=&totSpin, *p_energy_rnorm=&energy_rnorm;
    for(int i=0;i<ninterval;i++) totSpin+=initSpin[i];
    
    for(int i=0;i<nthermal*ninterval;i++) (*p_mcUpdate)(totOrbs, lattice, p_energy, p_totSpin); //thermalization

    double spin_i=0;
    double spin_j=0;
    double spin_ij=0;
    double totEnergy=0, totEnergy_rnorm=0;
    double E2=0, E2_rnorm=0;
    *p_energy=0;
    double M=0,M2=0,M4=0;
    double M_tmp=0,MdotM_tmp=0,M_tot=0;

    double spin_tot=0;
    
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
    for(int i=0;i<nsweep;i++){
        //printf("sweep%d/%d\n",i,nsweep);
        for(int j=0;j<ninterval;j++) (*p_mcUpdate)(totOrbs, lattice, p_energy, p_totSpin); // one sweep
        //printf("start statistics\n");
        //printf("spinFrame %d \n",iFrame);
        // record spin distribution
        if(spinFrame>0 & i%output_per_sweep==0){
            PyObject *spinDistribution=PyTuple_New(totOrbs);
            for(int j=0;j<totOrbs;j++) PyTuple_SetItem(spinDistribution, j, PyFloat_FromDouble(lattice[j].spin));
            PyTuple_SetItem(spinFrameData, iFrame, spinDistribution);
            iFrame+=1;
        }

        *p_energy_rnorm=0;
        for(int j=0;j<totOrbs;j++){ // calc. energy in renormalized system
            if(lattice[j].chosen>0) *p_energy_rnorm+=getCorrEnergy_rnorm(lattice+j)/2-(lattice+j)->h*(lattice+j)->spin;
        }
        //printf("calc. energy for renormalized lattice done\n");
        
        // spin statistics over space in each frame
        double spin_i_avg=0;
        double spin_j_avg=0;
        double corrAvg=0.0;
        for(int j=0;j<nLat;j++){
            double si_tmp=lattice[corrOrbPair[j*2+0]].spin;
            double sj_tmp=lattice[corrOrbPair[j*2+1]].spin;
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

        //if(h<1e-6){
            spin_i+=fabs(spin_i_avg)/nLat;//fabs(spin_i_avg)/nLat;
            spin_j+=fabs(spin_j_avg)/nLat;//fabs(spin_j_avg)/nLat;
        //}else{
        //    spin_i+=(spin_i_avg)/nLat;//fabs(spin_i_avg)/nLat;
        //    spin_j+=(spin_j_avg)/nLat;//fabs(spin_j_avg)/nLat;
        //}
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
    Data=PyTuple_New(11);
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
    PyTuple_SetItem(Data, 10, spinFrameData);
    //printf("before clib return data\n");
    return Data;
    
}

static PyMethodDef module_methods[] = {
    {"MCMainFunction", MCMainFunction, METH_VARARGS, "execute Monte Carlo sims. on Ising model"},
    {NULL, NULL, 0, NULL} // neccessary to tell python compiler stop here
};

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "isinglib",
    "Used to execute the Monte Carlo simulations of Ising model, that is, the O(1) spin model.",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_isinglib(void) {
    printf("Initializing isinglib...\n");
    PyObject* m;
    m= PyModule_Create(&moduledef);
    return m;
}