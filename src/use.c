// A couple of files that more than one thingy uses
#include<R.h>
#include<Rmath.h>
#include<Rinternals.h>
#include<stdlib.h>

#define eps 1.0e-32
#define EPSILON 0.0001   // Define your own tolerance
#define FLOAT_EQ(x,v) (((v - EPSILON) < x) && (x <( v + EPSILON)))

int diff(double *arr,int *len,double *retarr)
{
// this function takes a one-dimensional array arr and its length len, and returns the differenced
// vector retarr of length len-1

int i;
for(i=0;i<*len-1;i++)
{
	retarr[i] = arr[i+1]-arr[i];
}
return(0);

}

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

////////////////////// Some other functions //////////////////////////////////



double linearinterp(int n, double newx, double *a, double *b)
{
    double newvalue;
    int i;

//condition is where y lies between the two closests approximations to 
//it in the cal curve
   
    for(i=0; i<n-1; i++)
    {
        if (((newx >= a[i]) & (newx <= a[i+1])) | ((newx <= a[i]) & (newx >= a[i+1])))
        {
                newvalue = b[i] + ((newx-a[i])/(a[i+1]-a[i]))*(b[i+1]-b[i]);
                if(newx==a[i]) newvalue = b[i];
                return(newvalue);
                //break;
        }        
    }
  
  return(-999.0);
}

double Max2 (double a, double b)
{
	// find the max of 2 numbers
   double larger;
   if (a > b)
      larger = a;
   else
      larger = b;
   return larger;
}

//rtruncn function:
double rtruncn (double a, double b)
{
    double A, B, maxA, maxB, maxR, r2, r, th, u, x;
	int accept=0;
    
    A = atan(a);
    B = atan(b);
    
    maxA = exp(-pow(a,2)/4)/cos(A);
    maxB = exp(-pow(b,2)/4)/cos(B);
    maxR = Max2(maxA, maxB);

    if((a<1) && (b>-1)) maxR = exp(-0.25)*sqrt(2.0);

    while (accept!=1)
    {
        r2 = runif(0.0,1.0);
        r = sqrt(r2)*maxR;
        th = runif(A,B);
        u = r*cos(th);
        //v = r*sin(th);
        x = tan(th);
        accept = ((pow(x,2)) < (log(u)*-4));
    }
    return x;

}        

//truncated normal function:
double truncatedwalk (double old, double sd, double low, double high)
{
    double lowlimold, upplimold, y, newvalue;
    lowlimold = (low - old)/sd;
    upplimold = (high - old)/sd;
    y = rtruncn(lowlimold, upplimold);
    newvalue = old + sd*y;
           
    return newvalue;
}

//truncated normal ratio function:
double truncatedrat (double old, double sd, double low, double high, double newvalue)
{
    double lowlimold, upplimold, lowlimnew, upplimnew, plowold, puppold, plownew, puppnew, ratio;
    
    lowlimold = (low - old)/sd;
    upplimold = (high - old)/sd;
    lowlimnew = (low - newvalue)/sd;
    upplimnew = (high - newvalue)/sd;
    plowold = pnorm(lowlimold,0.0,1.0,1,0);
    puppold = pnorm(upplimold,0.0,1.0,1,0);
    plownew = pnorm(lowlimnew,0.0,1.0,1,0);
    puppnew = pnorm(upplimnew,0.0,1.0,1,0);
    ratio = (puppold - plowold)/(puppnew - plownew);
    return ratio;        
}

double Max(double *Numbers, int Count)
{
	// Find the maximum of a sequence of numbers
	double Maximum;
	Maximum = Numbers[0];

	for(int i = 0; i < Count; i++)
		if( Maximum < Numbers[i] )
			Maximum = Numbers[i];

	return Maximum;
}

double Min(double *Numbers, int Count)
{
	// Find the maximum of a sequence of numbers
	double Minimum;
	Minimum = Numbers[0];

	for(int i = 0; i < Count; i++)
		if( Minimum > Numbers[i] )
			Minimum = Numbers[i];

	return Minimum;
}

int seq(double from,double to,double len,double *sequence)
{
	// Create a sequence of numbers from 'from' to 'to' of length 'len'
	// Simple huh?
	
	double by = (to-from)/(len-1);
	int i;
	for(i=0;i<len;i++) 
		sequence[i] = from + i*by;

	return(0);

}


double UpdateMCMC(double logdiff,double newval,double oldval,double rat)
{
// Function to update MCMC when given log likelihoods

double u,mh;
u = runif(0.0,1.0);
mh = exp(logdiff)*rat;
if (u < mh)
{
	return(newval);
} else 
{
	return(oldval);
}

}

int fact(int number) {
	int temp;

	if(number <= 1) return 1;

	temp = number * fact(number - 1);
	return temp;
}

int sample(int *myseq, int nseq, double *pro) {
// take a sequence of values (in myseq) of length nseq, and output a single value 
// according to the weights in the prorow'th row of the proportions pro. 

int i,ans;
double u, *cumpro;
cumpro = (double *)calloc(nseq, sizeof(double));
if(cumpro==NULL) error("Can't allocate memory");

cumpro[0] = pro[0];
for(i=1;i<nseq;i++) cumpro[i] = pro[i]+cumpro[i-1];
u = runif(0.0,1.0);
// See where u lies in cumpro
ans=0;
for(i=0;i<nseq;i++) {
	if(u<cumpro[i]) {
		ans = i;
		break;
	}			
}

free(cumpro);
return(myseq[ans]);
}

void GaussJordan(int N, double **b, double **y)
{
//use generate_identity function to make y the indentity. Then call GaussJordan
//with a the matrix to invert. a gets turned into the identity, so make a copy
//if you need to keep it. N is the dimension of the matrix. 
	int             c, r, r_max, i, j;
	double          temp, v_v, v_max, factor = 0.0;

// Make it non-destructive
double **a;
a = (double **)calloc(N, sizeof(double *));
if(a==NULL) error("Can't allocate memory");
for (i=0;i<N;i++) a[i] = (double *)calloc(N, sizeof(double));

for(i=0;i<N;i++) {
  for(j=0;j<N;j++) {
    a[i][j] = b[i][j];
  }
}

	/* Loop over all columns of A */
	for (c = 0; c < N; c++)
	{			/* Find row with the maximum value absolute
				 * value. */
		r_max = c;
		v_max = fabs(a[c][c]);
		for (r = c + 1; r < N; r++)
		{
			v_v = fabs(a[r][c]);
			if (v_v > v_max)
			{
				r_max = r;
				v_max = v_v;
			}
		} /* Switch rows if necessary */ if (r_max != c)
		{
			for (j = c; j < N; j++)
			{
				temp = a[c][j];
				a[c][j] = a[r_max][j];
				a[r_max][j] = temp;
			} for (j = 0; j < N; j++)
			{
				temp = y[c][j];
				y[c][j] = y[r_max][j];
				y[r_max][j] = temp;
			}
		}		/* Rescale current row so that diagonal
		    element is 1 */ factor = a[c][c];
		if (fabs(factor) <= eps)
		{
			error("singular matrix\n");
		} for (j = c; j < N; j++)
		{
			a[c][j] /= factor;
		} for (j = 0; j < N; j++)
		{
			y[c][j] /= factor;
		} /* Subtract current row from all rows below. */
		for (r = c + 1; r < N; r++)
		{
			factor = a[r][c];
			for (j = c; j < N; j++)
				a[r][j] -= factor * a[c][j];
			for (j = 0; j < N; j++)
				y[r][j] -= factor * y[c][j];
		}
	
	}
	/*loop over all rows from 1 to end*/
	for (c = 1; c < N; c++)
		{
		/*Subtract current row from all rows above.*/
		for (j = 0; j < N-c; j++)
			{
			factor = a[j][N-c];
			for (r = 0; r < N; r++)
				{
				a[j][r] =a[j][r] - factor * a[N-c][r];
				y[j][r] =y[j][r] - factor * y[N-c][r];
				}
			}
		}
	for (i=0;i<N;i++) free(a[i]);
	free(a);
	return;
}

void cholesky(double **A, int n)
{
// Alternative cholesky decomposition from numerical recipes book
// NOTE: DESTRUCTIVE!

int i,j,k;
double sum;

for(i=0;i<n;i++) {
	for (j=i;j<n;j++) {
		sum=A[i][j];
		for (k=i-1;k>=0;k--) sum -= A[i][k]*A[j][k];
		if (i == j) {
			if (sum <= 0.0)	{
				error("Cholesky failed");
			}
			A[i][i]=sqrt(sum);
		} else A[j][i]=sum/A[i][i];
	}
}
// Set all uppder diagonals equal to zero
for (i=0;i<n;i++) for (j=0;j<i;j++) A[j][i] = 0.;

// Finally tranpose it - don't need this unless I'm doing random normal simulation
/*
double tmp;
for (i=0;i<n;i++) {
	for (j=0;j<=i;j++) {
		tmp = A[i][j];
		A[i][j] = A[j][i];
		A[j][i] = tmp;
	}
}
*/

}

double logdet(double **A, int n)
{
// Find log determinant using cholesky decomposition
// NON-DESTRUCTIVE

double sum = 0.0,**B;
int i,j;

B = (double **)calloc(n, sizeof(double *));
for (i=0;i<n;i++) B[i] = (double *)calloc(n, sizeof(double));
if(B==NULL) error("Can't allocate memory");

for (i=0;i<n;i++) {
	for (j=0;j<=i;j++) {
		B[i][j] = A[i][j];
		B[j][i] = A[j][i];
	}
}

cholesky(B,n);
for(i=0;i<n;i++) sum += log(B[i][i]);

for(i=0;i<n;i++) free(B[i]);
free(B);

return 2.0*sum;

}

double dlognormal(double *x, double *mean, double **var,int len)
{

// The log density of the multivariate normal distribution
// (-nrow(y)/2)*log(2*pi)-0.5*log(det(cov))-0.5*t(y-mean)%*%solve(cov)%*%(y-mean)
double dens,logdeter,**invvar,bigsum=0.0;
int i,j;

invvar = (double **)calloc(len, sizeof(double *));
for (i=0;i<len;i++) invvar[i] = (double *)calloc(len, sizeof(double));
if(invvar==NULL) error("Can't allocate memory");

for(i=0;i<len;i++) invvar[i][i] = 1.0;
for(i=0;i<len;i++) {
	for(j=0;j<len;j++) {
		invvar[i][j] = 1.0;
		invvar[j][i] = 1.0;
	}
}

logdeter = logdet(var,len);
GaussJordan(len,var,invvar);

for(i=0;i<len;i++) {
  for(j=0;j<len;j++) {
    bigsum += (x[j]-mean[j])*(x[i]-mean[i])*invvar[i][j];
  }
}

if(R_FINITE(logdeter)==FALSE) {
	dens = -1/eps;
	error("Determinant non-finite \n");
} else {
	dens = -(((double)len)/2) * log(2*PI) - 0.5*logdeter - 0.5*bigsum;
}

for(i=0;i<len;i++) free(invvar[i]);
free(invvar);
return(dens);

}


void ranmvn(int N, double *mu, double **P, double *v) {
/*
 You want a random vector from N(mu, S).  P is the transpose
 of the cholesky decomposition of S (if S = PP', you send me P').
 mu is N x 1
 P is N x N
 v is the vector you want to output
 
 This function will fillup the 1st N elements of v so
 as to be a single random draw from the MVN(mu, S) density.
 
 Hard limit: N cannot be > 1000.
 */
	if(N>1000) error("Climate dimension too large! \n");
	int i, j; 
	double sum, *tmp;
	tmp = (double *)calloc(1000, sizeof(double));
	if(tmp==NULL) error("Can't allocate memory");

	for (i=0;i<N;i++) tmp[i] = rnorm(0,1);
	for (i=0;i<N;i++) {
		sum = 0.0;
		for(j=0;j<N;j++) sum += tmp[j]*(P[i][j]);
		v[i] = sum+mu[i];
	}
	free(tmp);
}

// My function starts here, take a mean vector and variance matrix and produce a single random sample
void rmvnorm(double *out, double *mean, double **var,int length) {

	// V is the variance matrix, transchol will be the transpose of the cholesky decomposition, 
	int i,j;
	double **transchol;

	transchol = (double **)calloc(length, sizeof(double *));
	for (i=0;i<length;i++) transchol[i] = (double *)calloc(length, sizeof(double));
	if(transchol==NULL) error("Can't allocate memory");

	for(i=0;i<length;i++) {
		for(j=0;j<=i;j++) {
			transchol[i][j] = var[i][j];
			transchol[j][i] = var[i][j];
		}
	}

	cholesky(transchol,length);
	
	//for(i=0;i<length;i++) mu[i] = mean[i];
	ranmvn(length, mean, transchol, out);
	for (i=0;i<length;i++) free(transchol[i]);
	free(transchol);
}

void CreateTplusW(double *RWprecs, double *vars, int dim, double **out) {

int i;
// Set the diagonals 
for(i=0;i<dim;i++) out[i][i] = 1/vars[i];
for(i=0;i<dim-1;i++) out[i][i] += RWprecs[i];
for(i=1;i<dim;i++) out[i][i] += RWprecs[i-1];

// Set the off diagonals
for(i=0;i<dim-1;i++) out[i][i+1] = -RWprecs[i];
for(i=0;i<dim-1;i++) out[i+1][i] = -RWprecs[i];

}



double dlinvgauss(double x, double mu, double lambda) {
// Log inverse gaussian density
double dens;
// Note: this is the R version of the IG distribution (different to the one used in the NIG paper)
dens = (log(lambda) - log(2 * PI) - 3 * log(x))/2 - lambda * pow(x - mu,2)/(2 * pow(mu,2) * x);
return(dens);
}

void trisolve (int n, double *a, double *b, double *c, double *v, double *x)
{
    /**
     * n - number of equations
     * a - sub-diagonal (means it is the diagonal below the main diagonal) -- indexed from 1..n-1
     * b - the main diagonal
     * c - sup-diagonal (means it is the diagonal above the main diagonal) -- indexed from 0..n-2
     * v - right part
     * x - the answer
     Ammended to make it non-destructive of b and v
     */
    int i;
    double *b2,*v2;
    b2 = (double *)calloc(n, sizeof(double));
    v2 = (double *)calloc(n, sizeof(double));
    if(b2==NULL) error("Can't allocate memory");
    if(v2==NULL) error("Can't allocate memory");
    for(i=0;i<n;i++) {
        b2[i] = b[i];
        v2[i] = v[i];
    }
    
    for (int i = 1; i < n; i++)
    {
        double m = a[i]/b2[i-1];
        b2[i] = b2[i] - m*c[i-1];
        v2[i] = v2[i] - m*v2[i-1];
    }
    
    x[n-1] = v2[n-1]/b2[n-1];
    
    for (int i = n - 2; i >= 0; i--)
        x[i]=(v2[i]-c[i]*x[i+1])/b2[i];
        
    free(b2);
    free(v2);
}

void maketri(double *v, int n, double *D, double*upp,double *dia, double *low) {

    // Function to make tri-diagonal parts of Q/V matrix
    int l;
    upp[n-1] = 0.0;
    for(l=0;l<(n-1);l++) upp[l] = -1/(v[l]);
    dia[0] = 1/(v[0])+D[0];
    dia[n-1] = 1/(v[n-2])+D[n-1];
    for(l=1;l<(n-1);l++) dia[l] = 1/(v[l-1])+1/(v[l])+D[l];
    low[0] = 0.0;
    for(l=1;l<n;l++) low[l] = -1/(v[l-1]);    

}

int samplegrid(int gridsize,int oldindex,int step) {
    // Function to move from old value to new value along a grid where the size of the move is controlled by the variable step
    int lowindex,highindex,ans;
    double U;
    
    lowindex = 0;
    if(oldindex-step>0) lowindex = oldindex-step;
    highindex = gridsize-1;
    if(oldindex+step<highindex) highindex= oldindex+step;
    U = runif(0,1);
    ans = (int)ftrunc((1.0+(double)highindex-(double)lowindex)*U+(double)lowindex);
    
    return(ans);
    
}

double dlinvgauss2(double x, double mu, double phi) {
// Alternative version of log inverse gaussian density
// Note: this is the Karlis (IG2) version of the IG distribution (different to the one used in R)
double dens;
dens = 0.5*log(mu)+0.5*log(phi)-1.5*log(x)+phi-0.5*phi*(x/mu+mu/x)-0.5*log(2*PI);
return(dens);
}

void CreateUs(int n, double *vj, double phipar2, double u1, double u2, double u3) {

// Function to create u-values used by main function for updating mu and phi
    double meanv=0, meanvr=0;
    int i;
    for(i=0;i<(n-1);i++) {
        meanv += vj[i];
        meanvr += 1/vj[i];
    }
    meanv = meanv/(n-1);
    meanvr = meanvr/(n-1);
        
    u1 = ((double)n-1)*meanv;
    u2 = ((double)n-1)-phipar2;
    u3 = ((double)n-1)*meanvr;

    //Rprintf("u1=%lf, u2=%lf u3=%lf \n", u1,u2,u3);

}

void CholTriDiag(double *alpha, double *beta, int n, double *delta, double *l) {
// Function to create cholesky decomposition of tri-diagonal matrix
// alpha and beta are the diagonals and off-diagonals of the input matrix respectively
// delta and l are the diagonals and off-diagonals of the output matrix respectively
    int i;
    delta[0] = alpha[0];
    l[0] = beta[0]/sqrt(delta[0]);
    for(i=1;i<n-1;i++) {
        delta[i] = alpha[i]-pow(beta[i-1],2)/delta[i-1];
        l[i] = beta[i]/sqrt(delta[i]);
    }
    delta[n-1] = alpha[n-1]-pow(beta[n-2],2)/delta[n-2];
    for(i=0;i<n;i++) delta[i] = sqrt(delta[i]);
}

double logdetTriDiag(double *diag, double *upper, double *lower, int n) {
// Function to calculate the determinant for a tri-diagonal matrix
// Idea taken from http://en.wikipedia.org/wiki/Tridiagonal_matrix
// diag is the diagonal, upper is the upper diagonal, lower is the lower
    int i;
    double *f,logdet;
    f = (double *)calloc(n, sizeof(double));
    if(f==NULL) error("Can't allocate memory");
    f[0] = 1;
    f[1] = diag[0];
    for(i=2;i<n;i++) {
        f[i] = diag[i-1]*f[i-1]-lower[i-2]*upper[i-2]*f[i-2];
    }
    logdet = log(diag[n-1]*f[n-1]-lower[n-2]*upper[n-2]*f[n-2]);
    free(f);
    return(logdet);
}


    