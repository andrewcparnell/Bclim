// This function runs the Bclim MCMC functions

#include<R.h>
#include<Rmath.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include"use.h"
#include <R_ext/Utils.h>

#define LARGE 1.0e20
   // Define your own tolerance

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
////////////////////// Main MCMC function //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

void BclimMCMC3D(int *G,int *n, int *m,int *nchrons, double *MixPro, double *MixMeans, double *MixPrec, char **ChronsFILE, int *iterations, int *burnin, int *thinby, int *reportevery, double *vmhsd, double *vstart, int *Zstart, double *phi1start, double *phi2start, double *vstore, int *zstore, double *chronstore, double *cstore, double *phi1store, double *phi2store, double *phi1dlmean, double *phi1dlsd, double *phi2dlmean, double *phi2dlsd, double *phi1mhsd, double *phi2mhsd)
{

// G is number of mixture groups, n is number of layers, m is number of climate dimensions (always 3 here), nchrons is number of chronologies
// MixPro is mixture proportions, MixMeans is mixture means, MixPrec is mixture precisions
// ChronsFILE is the path to the chronologies file
// iterations is number of iterations, burnin is size of burnin, thinby is thinning amount, reportevery is how often to report
// vmhsd is the standard deviation of the truncated random walk for proposing new values of v
// vstart and Zstart are starting values for v and Z respectively
// phi1 and phi2 are 3-vectors of the values used for phi1 and phi2
// vstore is variance output, chronstore are used chronologies, cstore is used for climates

// Declare indicators
int i,j,k,l,j2,j3,iter;
	
// Get a new seed
GetRNGstate();

// Set up arrays to enable simpler mixture read-ins
double ***MyMixMean,**MyMixPrec,**MyMixPro,*CurrMixPrec;
MyMixMean = (double ***)calloc(*n, sizeof(double **));
for(i=0;i<*n;i++) MyMixMean[i] = (double **)calloc(*m,sizeof(double *));
for(i=0;i<*n;i++) for(j=0;j<*m;j++) MyMixMean[i][j] = (double *)calloc(*G, sizeof(double));

MyMixPrec = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) MyMixPrec[i] = (double *)calloc(*G, sizeof(double));
MyMixPro = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) MyMixPro[i] = (double *)calloc(*G, sizeof(double));
CurrMixPrec = (double *)calloc(*n, sizeof(double));
if(MyMixMean==NULL) error("Can't allocate memory");
if(MyMixPrec==NULL) error("Can't allocate memory");
if(MyMixPro==NULL) error("Can't allocate memory");
if(CurrMixPrec==NULL) error("Can't allocate memory");

// Turn MixPrec, MixMeans and MixPro into arrays for easy lookup
for(i=0;i<*n;i++) for(k=0;k<*G;k++) MyMixPrec[i][k] = MixPrec[k*(*n)+i];
//for(i=0;i<*n;i++) for(k=0;k<*G;k++) Rprintf("%lf \n",MyMixPrec[i][k]);
for(i=0;i<*n;i++) {
	for(j=0;j<*m;j++) {
		for(k=0;k<*G;k++) {
			MyMixMean[i][j][k] = MixMeans[k*(*m)*(*n)+j*(*n)+i];
      //Rprintf("%lf \n",MyMixMean[i][j][k]);
		}
	}
}
for(i=0;i<*n;i++) for(j=0;j<*G;j++) MyMixPro[i][j] = MixPro[j*(*n)+i];

// Set up matrices
double **v,*phi1,*phi2,*mixpro,*choldiag,*chollower,*rannormal,*cholzeros,*mvnvar;
int *Z,*Zstar,*mixseq;
v = (double **)calloc(*n-1, sizeof(double *));
for(i=0;i<*n-1;i++) v[i] = (double *)calloc(*m, sizeof(double));
phi1 = (double *)calloc(*m, sizeof(double));
phi2 = (double *)calloc(*m, sizeof(double));
Z = (int *)calloc(*n, sizeof(int));
Zstar = (int *)calloc(*n, sizeof(int));
mixpro = (double *)calloc(*G, sizeof(double));
for(i=0;i<*G;i++) mixpro[i]=1/((double)*G);
mixseq = (int *)calloc(*G, sizeof(int));
for(i=0;i<*G;i++) mixseq[i]=i;
choldiag = (double *)calloc(*n, sizeof(double));
chollower = (double *)calloc(*n, sizeof(double));
rannormal = (double *)calloc(*n, sizeof(double));
cholzeros = (double *)calloc(*n, sizeof(double));
mvnvar = (double *)calloc(*n, sizeof(double));

if(v==NULL) error("Can't allocate memory");
if(phi1==NULL) error("Can't allocate memory");
if(phi2==NULL) error("Can't allocate memory");
if(Z==NULL) error("Can't allocate memory");
if(Zstar==NULL) error("Can't allocate memory");
if(mixseq==NULL) error("Can't allocate memory");
if(choldiag==NULL) error("Can't allocate memory");
if(chollower==NULL) error("Can't allocate memory");
if(rannormal==NULL) error("Can't allocate memory");
if(cholzeros==NULL) error("Can't allocate memory");
if(mvnvar==NULL) error("Can't allocate memory");

// Give starting values
for(i=0;i<*n-1;i++) for(j=0;j<*m;j++) v[i][j] =  vstart[i];
for(i=0;i<*m;i++) phi1[i] = phi1start[i];
for(i=0;i<*m;i++) phi2[i] = phi2start[i];
for(i=0;i<*n;i++) Z[i] = Zstart[i]-1;

// Set up the matrices I'll need;
double **M, **R, **DM, **Mstar, **Rstar, **DMstar;
M = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) M[i] = (double *)calloc(*m, sizeof(double)); // n rows, m columns - hopefully!
R = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) R[i] = (double *)calloc(*m, sizeof(double));
DM = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) DM[i] = (double *)calloc(*m, sizeof(double));
Mstar = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) Mstar[i] = (double *)calloc(*m, sizeof(double));
Rstar = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) Rstar[i] = (double *)calloc(*m, sizeof(double));
DMstar = (double **)calloc(*n, sizeof(double *));
for(i=0;i<*n;i++) DMstar[i] = (double *)calloc(*m, sizeof(double));

// Some final vectors I need
double *MDR,*MDM,*MDRstar,*MDMstar,*D,*Wupper,*Wdiag;
D = (double *)calloc(*n, sizeof(double));
MDR = (double *)calloc(*m, sizeof(double));
MDM = (double *)calloc(*m, sizeof(double));
MDRstar = (double *)calloc(*m, sizeof(double));
MDMstar = (double *)calloc(*m, sizeof(double));

if(M==NULL) error("Can't allocate memory");
if(D==NULL) error("Can't allocate memory");
if(R==NULL) error("Can't allocate memory");
if(DM==NULL) error("Can't allocate memory");
if(Mstar==NULL) error("Can't allocate memory");
if(Rstar==NULL) error("Can't allocate memory");
if(DMstar==NULL) error("Can't allocate memory");

// Set up somewhere to store output
double **allchrons;
allchrons = (double **)calloc(*nchrons, sizeof(double *));
for(i=0;i<*nchrons;i++) allchrons[i] = (double *)calloc(*n, sizeof(double));
if(allchrons==NULL) error("Can't allocate memory");

// Read chronologies into matrix to make it easier to randomly select one
FILE *chrons;
chrons = fopen(*ChronsFILE,"r");
if(chrons==NULL) error("Error: can't open chronologies file %s.\n",*ChronsFILE);		
double tmp;	
for(i=0;i<*nchrons;i++) for(j=0;j<*n;j++) tmp=fscanf(chrons,"%lf",&allchrons[i][j]);                       
fclose(chrons);

// Print one line of chronologies if required
//for(i=1;i<*n;i++) Rprintf("%lf \n",allchrons[0][i]);

// Create chron differences
double *diffchron,*currchron;
diffchron = (double *)calloc(*n-1, sizeof(double));
currchron = (double *)calloc(*n, sizeof(double));
if(diffchron==NULL) error("Can't allocate memory");
if(currchron==NULL) error("Can't allocate memory");

// Proposal vector for Z;
double *eqprob;
eqprob = (double *)calloc(*G, sizeof(double));
for(k=0;k<*G;k++) eqprob[k] = (double)1 / *G;

// Set up MCMC variables I'll need
double *tempv,*tempDM,*tempR,*tempDMstar,*tempRstar,*tempvstar,*triupper,*tridiag,*trilower,*logdetratio,vstar,vstarrat,logRv,logRz,*tempc,phi1star,phi2star,phi1starrat,phi2starrat,logRphi1,logRphi2,*B,*S,tBS;
tempv = (double *)calloc(*n-1, sizeof(double));
tempvstar = (double *)calloc(*n-1, sizeof(double));
tempDM = (double *)calloc(*n, sizeof(double));
tempR = (double *)calloc(*n, sizeof(double));
tempDMstar = (double *)calloc(*n, sizeof(double));
tempRstar = (double *)calloc(*n, sizeof(double));
logdetratio = (double *)calloc(*m, sizeof(double));
triupper = (double *)calloc(*n, sizeof(double));
tridiag = (double *)calloc(*n, sizeof(double));
trilower = (double *)calloc(*n, sizeof(double));
tempc = (double *)calloc(*n, sizeof(double));
B = (double *)calloc(*n, sizeof(double));
S = (double *)calloc(*n, sizeof(double));
int accept,countchron=0,countv=0,countz=0,currentpos,countc=0,countphi1=0,countphi2=0;

// Start off the chronologies
for(i=0;i<*n;i++) currchron[i] = allchrons[0][i];
diff(currchron,n,diffchron);

// Create the matrices I need later
for(j=0;j<*m;j++) {
    // Create M,D and DM and MDM, etc
    MDM[j] = 0;
    MDMstar[j] = 0;
    for(i=0;i<*n;i++) {
        M[i][j] =  MyMixMean[i][j][Z[i]];
        D[i] = MyMixPrec[i][Z[i]];
        DM[i][j] = D[i]*M[i][j];
        Mstar[i][j] =  MyMixMean[i][j][Z[i]];
        DMstar[i][j] = D[i]*M[i][j];
        MDM[j] += M[i][j]*DM[i][j];
        MDMstar[j] += M[i][j]*DM[i][j];
     }
}

// Start iterations loop
double progress=0;
for(iter=0;iter<*iterations;iter++) {


	// Check if someone's pressed escape
	R_CheckUserInterrupt();

    // Report progress
	if(iter%*reportevery==0) {
	    progress = (double) 100*iter/ *iterations;
        Rprintf("\r");
        Rprintf("Optimising parameters: %3.2f %% completed",progress);
        //Rprintf("Completed: %i ",iter);
	    Rprintf("\r");
	    R_FlushConsole();
	}
    
    // Store stuff
    if((iter%*thinby==0) & (iter>=*burnin)) {
        
        currentpos = (int)((iter-*burnin)/(*thinby));
        // Sample a chronology
        for(i=0;i<*n;i++) currchron[i] = allchrons[(int)currentpos][i];
        diff(currchron,n,diffchron);
		
        // Store the used chronology
        for(i=0;i<*n;i++) chronstore[countchron+i] = currchron[i];
        countchron += *n;

        // Store current z values
		for(i=0;i<*n;i++) zstore[countz+i] = Z[i];
		countz += *n;
        
        // Store phi1 and phi2
        for(j=0;j<*m;j++) {
            phi1store[countphi1+j] = phi1[j];
            phi2store[countphi2+j] = phi2[j];
        }
        countphi1 += *m;
        countphi2 += *m;
        
        for(j=0;j<*m;j++) {
            
            for(i=0;i<*n-1;i++) vstore[countv+i] = v[i][j];
            countv += *n-1;

            // Generate n independent random normals
            for(l=0;l<*n;l++) rannormal[l] = rnorm(0,1);

            for(l=0;l<*n-1;l++) tempv[l] = v[l][j];
            for(l=0;l<*n;l++) tempDM[l] = DM[l][j];
            for(l=0;l<*n;l++) tempR[l] = R[l][j];
            
            // Create the tri-diagonal martix Q, stored in triupper, tridiag and trilower
            maketri(tempv,*n,D,triupper,tridiag,trilower);

            // Create cholesky decomposition of the tridiag
            CholTriDiag(tridiag,triupper,*n,choldiag,chollower);

            // Solve the cholesky decomposition
            trisolve(*n, chollower,choldiag, cholzeros, rannormal, mvnvar);

            // Get the mean
            trisolve(*n, trilower,tridiag, triupper, tempDM, tempR);

            // Output the climates
            for(i=0;i<*n;i++) tempc[i] = tempR[i] + mvnvar[i];

            // Store the climates in the right place
            for(i=0;i<*n;i++) cstore[countc+i] = tempc[i];
            countc += *n;
            // Need to be careful with c as now we're storing them in a weird order
            
        // End of j loop through climate dimensions
        }
        
    
    // End of storage loop
    }

    // Loop through climate dimensions to update v
    for(j=0;j<*m;j++) {
            
        for(l=0;l<*n-1;l++) tempv[l] = v[l][j];
        for(l=0;l<*n;l++) tempDM[l] = DM[l][j];
        for(l=0;l<*n;l++) tempR[l] = R[l][j];
        
        // Create the tri-diagonal martix Q, stored in triupper, tridiag and trilower
        maketri(tempv,*n,D,triupper,tridiag,trilower);
        // Find R from solution given in appendix to paper
        trisolve(*n,trilower,tridiag,triupper,tempDM,tempR);
        // Finally create MDR
        MDR[j] = 0;
        for(l=0;l<*n;l++) MDR[j] += tempDM[l]*tempR[l];
        
        // Sample a new v
        for(i=0;i<*n-1;i++) {
            
            // For some reason I found it easier to accept new values on the volatility scale
            vstar = truncatedwalk(v[i][j], *vmhsd, 0, LARGE);
            vstarrat = truncatedrat(v[i][j], *vmhsd, 0, LARGE, vstar);
            for(l=0;l<*n-1;l++) tempvstar[l] = tempv[l];
            tempvstar[i] = vstar;
            
            // Find log ratio of determinant
            for(l=0;l<*n;l++) B[l] = 0.0;
            B[i] = 1.0;
            B[i+1] = -1.0;
            trisolve (*n, trilower,tridiag, triupper, B, S);
            tBS = S[i]-S[i+1];
            logdetratio[j] = log(1+(1/vstar-1/v[i][j])*tBS);
            
            // Calculate Rstar
            maketri(tempvstar,*n,D,triupper,tridiag,trilower);
            trisolve(*n,trilower,tridiag,triupper,tempDM,tempRstar);

            //Calculate MDRstar
            MDRstar[j] = 0;
            for(l=0;l<*n;l++) MDRstar[j] += tempDM[l]*tempRstar[l];

            logRv = -0.5*log(vstar/tempv[i]) -0.5*logdetratio[j] + 0.5*MDRstar[j] - 0.5*MDR[j] + dlinvgauss2(vstar,phi1[j]*diffchron[i],phi2[j]*diffchron[i]) - dlinvgauss2(tempv[i],phi1[j]*diffchron[i],phi2[j]*diffchron[i]);
            
            //Rprintf("logRv=%lf \n",logRv);
            
            //Rprintf("i=%i, j=%i, v=%lf, vstar=%lf, vstarrat=%lf logRv=%lf logdetratio[j]=%lf MDRstar[j]=%lf MDR[j]=%lf dligstar=%lf dlig=%lf\n",i,j,v[i][j],vstar,vstarrat,logRv,logdetratio[j],MDRstar[j],MDR[j],dlinvgauss2(vstar,phi1[j]*diffchron[i],phi2[j]*diffchron[i]),dlinvgauss2(tempv[i],phi1[j]*diffchron[i],phi2[j]*diffchron[i]));
            
            accept = (int)UpdateMCMC(logRv,1,0,vstarrat);
            if(accept==1) {
                v[i][j] = vstar;
                tempv[i] = vstar;
                maketri(tempv,*n,D,triupper,tridiag,trilower);
                for(l=0;l<*n;l++) tempR[l] = tempRstar[l];
                for(l=0;l<*n;l++) R[l][j] = tempRstar[l];
                MDR[j] = MDRstar[j];
            }
            
        } // End of i loop for v update

    // End of loop through climate dimensions
    }
    
    
    // Sample a new z
    for(i=0;i<*n;i++) {
        for(l=0;l<*n;l++) Zstar[l] = Z[l];
        Zstar[i] = sample(mixseq, *G, eqprob);
        
        for(j=0;j<*m;j++) {
            for(l=0;l<*n;l++) Mstar[l][j] =  MyMixMean[l][j][Zstar[l]];
            for(l=0;l<*n;l++) DMstar[l][j] = D[l]*Mstar[l][j];
            MDMstar[j] = MDM[j] - M[i][j]*DM[i][j] + Mstar[i][j]*DMstar[i][j];
            
            // Create Rstar
            for(l=0;l<*n-1;l++) tempv[l] = v[l][j];
            for(l=0;l<*n;l++) tempDM[l] = DMstar[l][j];
            for(l=0;l<*n;l++) tempR[l] = R[l][j];

            maketri(tempv,*n,D,triupper,tridiag,trilower);
            trisolve(*n,trilower,tridiag,triupper,tempDM,tempR);
           
            MDRstar[j] = 0;
            for(l=0;l<*n;l++) MDRstar[j] += tempDM[l]*tempR[l];            
        }
                
        logRz = log(MyMixPro[i][Zstar[i]]/MyMixPro[i][Z[i]]);
        for(j=0;j<*m;j++) {
            logRz += - 0.5*MDMstar[j] + 0.5*MDRstar[j] + 0.5*MDM[j] - 0.5*MDR[j];
        }
        
        accept = (int)UpdateMCMC(logRz,1,0,1);
        if(accept==1) {
            Z[i] = Zstar[i];
            for(j=0;j<*m;j++) {
                for(l=0;l<*n;l++) M[l][j] = Mstar[l][j];
                for(l=0;l<*n;l++) DM[l][j] = DMstar[l][j];
                MDM[j] = MDMstar[j];
                MDR[j] = MDRstar[j];
            
                for(l=0;l<*n-1;l++) tempv[l] = v[l][j];
                for(l=0;l<*n;l++) tempDM[l] = DMstar[l][j];
                for(l=0;l<*n;l++) tempR[l] = R[l][j];
                maketri(tempv,*n,D,triupper,tridiag,trilower);
                trisolve(*n,trilower,tridiag,triupper,tempDM,tempR);
                for(l=0;l<*n;l++) R[l][j] = tempR[l];
            }
        }
    
    // End of loop for z
	}
    
    // Sample new phi1 and phi2
    for(j=0;j<*m;j++) {
        phi1star = truncatedwalk(phi1[j], *phi1mhsd, 0, LARGE);
        phi1starrat = truncatedrat(phi1[j], *phi1mhsd, 0, LARGE, phi1star);

        logRphi1 = dlnorm(phi1star,phi1dlmean[j],phi1dlsd[j],1) - dlnorm(phi1[j],phi1dlmean[j],phi1dlsd[j],1);
        for(i=0;i<*n-1;i++) logRphi1 += dlinvgauss2(v[i][j],phi1star*diffchron[i],phi2[j]*diffchron[i]) - dlinvgauss2(v[i][j],phi1[j]*diffchron[i],phi2[j]*diffchron[i]);
        //for(i=0;i<*n-1;i++) Rprintf("i=%i, v=%lf, first=%lf, second=%lf diffchron=%lf \n",i,v[i][j],dlinvgauss2(v[i][j],phi1star*diffchron[i],phi2[j]*diffchron[i]), dlinvgauss2(v[i][j],phi1[j]*diffchron[i],phi2[j]*diffchron[i]),diffchron[i]);
        
        //Rprintf("logRphi1=%lf j=%i, phi1star=%lf, phi1starrat=%lf, phi1[j]=%lf, dlnorm(phi1star,phi1dlmean[j],phi1dlsd[j],1)=%lf, dlnorm(phi1[j],phi1dlmean[j],phi1dlsd[j],1)=%lf, \n",logRphi1,j,phi1star,phi1starrat,phi1[j],dlnorm(phi1star,phi1dlmean[j],phi1dlsd[j],1),dlnorm(phi1[j],phi1dlmean[j],phi1dlsd[j],1));
        phi1[j] = UpdateMCMC(logRphi1,phi1star,phi1[j],phi1starrat);

        phi2star = truncatedwalk(phi2[j], *phi2mhsd, 0, LARGE);
        phi2starrat = truncatedrat(phi2[j], *phi2mhsd, 0, LARGE, phi2star);

        logRphi2 = dlnorm(phi2star,phi2dlmean[j],phi2dlsd[j],1) - dlnorm(phi2[j],phi2dlmean[j],phi2dlsd[j],1);
        for(i=0;i<*n-1;i++) logRphi2 += dlinvgauss2(v[i][j],phi1[j]*diffchron[i],phi2star*diffchron[i]) - dlinvgauss2(v[i][j],phi1[j]*diffchron[i],phi2[j]*diffchron[i]);

        //Rprintf("logRphi2=%lf \n",logRphi2);
        phi2[j] = UpdateMCMC(logRphi2,phi2star,phi2[j],phi2starrat);
    
    
    // End of j loop updating phi1 and phi2
    }
    
    
// End of iterations loop
}

// Change the RNG state
PutRNGstate();

Rprintf("\r");
R_FlushConsole();
Rprintf("Completed: 100.00 %%");
Rprintf("\n");
R_FlushConsole();

// End of function
}
