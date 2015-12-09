// Header file for use.c

int diff(double *arr,int *len,double *retarr);
int compare (const void * a, const void * b);
double linearinterp(int n, double newx, double *a, double *b);
double Max2 (double a, double b);
double rtruncn (double a, double b);
double truncatedwalk (double old, double sd, double low, double high);
double truncatedrat (double old, double sd, double low, double high, double newvalue);
double Max(double *Numbers, int Count);
double Min(double *Numbers, int Count);
int seq(double from,double to,double len,double *sequence);
double UpdateMCMC(double logdiff,double newval,double oldval,double rat);
int fact(int number);
int sample(int *myseq, int nseq, double *pro);
void GaussJordan(int N, double **a, double **y);
void ranmvn(int N, double *mu, double **P, double *v);
void rmvnorm(double *out, double *mean, double **var,int length);
void cholesky(double **A, int n); // DESTRUCTIVE
double logdet(double **A, int n); // NON-DESTRUCTIVE
void CreateTplusW(double *RWprecs, double *precs, int dim, double **out);
double dlinvgauss(double x, double mu, double lambda);
void trisolve (int n, double *a, double *b, double *c, double *v, double *);
void maketri(double *v, int n, double *D, double*upp,double *dia, double *low);
int samplegrid(int gridsize,int oldindex,int step);
double dlinvgauss2(double x, double mu, double phi);
void CreateUs(int n, double *vj, double phipar2, double u1, double u2, double u3);
void CholTriDiag(double *alpha, double *beta, int n, double *delta, double *l);
double logdetTriDiag(double *a, double *b, double*c, int n);