/*------------------------------------------------------------------------------
              James Sweeney          Trinity College Dublin          sweeneja@tcd.ie

                                     20/7/2011

      -       Program to return a probability distribution on climate given a matrix of pollen counts
      -       Integration method used is laguerre quadrature - (20 pts?)
      -       3D climate used
-------------------------------------------------------------------------------*/

#include <Rmath.h>
#include <Rinternals.h>
#include <R.h>

/*------------------------------------------------------------------------------*/
//function for evaluating ZINB likelihood - returns log values

double ZINB(double y, double x,double alpha, double delta)
{
    double q;double ans; double mu;
    q = pow((exp(x)/(1+exp(x))),alpha);
    mu=exp(x);

    if( y == 0 ){
        ans = log(1-q+ q*dnbinom_mu(y, delta, mu,0));
    }else{
        ans = log(q) + dnbinom_mu(y, delta, mu,1);
    }
    return ans;
}

void PalaeoRecon3D(int *nslices,double *Resultslices,int *pollen,double *MeanRespSurf,double *VarRespSurf, double *alpha,double *delta,double *Buffer, int *nquadpts, double *quadpts, double *quadprobs)
{

    /*--------------------------------------------------------------------------
    // nslices = #slices; dim(Resultslices) = nslices*2500
    // dim(MeanRespSurf) = dim(VarRespSurf) = 28*2500
    // beta = vector of regression parameters; predictor = predictor variables
    // alpha = ZI parameter; delta=overdispersion parameter
    // Buffer = locations where climate is feasible
    // quadpts = laguerre quadrature pts; quadprobs = laguerre quad weights
    --------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------*/
//definition of necessary parameters

                   int i=0,j=0,k=0,l=0;double sum=0,q_t=0,a=0,b=0;

/*------------------------------------------------------------------------------*/
//for loop where the fun starts

double progress=0;
for(i=0;i<*nslices;i++){                                                          //iterates through slices

    progress = (double) 100*i/ *nslices;
    Rprintf("\r");
    Rprintf("%3.2f%% completed",progress);
	Rprintf("\r");
	R_FlushConsole();
	//fflush(stdout);

    //Rprintf("\r %i \n",*nslices-i,"\r");
    //R_FlushConsole();
    // Check if someone's pressed escape
    R_CheckUserInterrupt();

       for(k=0;k<175616;k++){                                                       //iterates through elements of response surface i

        if(Buffer[k]==1){

           for (j=0;j<28;j++){                                                      //iterates through response surfaces
                 sum=0;                                                             // required for summation

                     a= MeanRespSurf[k+(j*175616)]-5*pow(VarRespSurf[k+(j*175616)],.5); b=MeanRespSurf[k+(j*175616)]+5*pow(VarRespSurf[k+(j*175616)],.5);

                  for(l=0;l<*nquadpts;l++){                                                     // iterates through quadrature points
                                  q_t= ((b-a)/2)*quadpts[l] + ((a+b)/2);
                                  sum=sum+exp(ZINB(pollen[i+j*nslices[0]],(q_t),alpha[j],delta[j])+dnorm(q_t,MeanRespSurf[k+(j*175616)],pow(VarRespSurf[k+(j*175616)],.5),1)+log((b-a)/2)+log(quadprobs[l]));                     // probabilites
                                  }

                           Resultslices[i*175616 + k]= Resultslices[i*175616 + k] + log(sum);
                       }
                    }
             }
}

    Rprintf("\r");
    R_FlushConsole();
    Rprintf("Completed: 100.00 %%");
    Rprintf("\n");
    R_FlushConsole();

/*------------------------------------------------------------------------------*/
//end of program
}
