#include "math.h"
#include "mex.h"   //--This one is required

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double  *x;//input
	double *y;//output
	
     //mwSize mrows,ncols;
    unsigned int mrows,ncols;
	double *tk;//tranzition levels
	
	unsigned int N_x,N_tk;
	
	unsigned int ii,jj;
	
	if (nrhs != 2) 
	{ 
		mexErrMsgTxt("Two input arguments required."); 
	} 
	else if (nlhs > 1) 
	{
		mexErrMsgTxt("Too many output arguments."); 
	} 
	else if (!mxIsNumeric(prhs[0])) 
	{
		mexErrMsgTxt("Argument must be numeric.");
	} 
	else if (!mxIsNumeric(prhs[1])) 
	{
		mexErrMsgTxt("Argument must be numeric.");
	} 
    
	if( ! ( (mxGetN(prhs[0]) == 1) || (mxGetM(prhs[0]) == 1)))
	{
		mexErrMsgTxt("Tranzition levels must be a vector.");
	}

	if( ! ( (mxGetN(prhs[1]) == 1) || (mxGetM(prhs[1]) == 1)))
	{
		mexErrMsgTxt("Input data must be a vector.");
	}
	
	if(mxGetN(prhs[1]) == 1) 
	{
		N_x = mxGetM(prhs[1]);
	}
	else
	{
		N_x = mxGetN(prhs[1]);
	}
	
//	printf("N_x = %d\n",N_x);
    
    mrows = mxGetM(prhs[1]);
    ncols =  mxGetN(prhs[1]);
        
	if(mxGetN(prhs[0]) == 1) 
	{
		N_tk = mxGetM(prhs[0]);
	}
	else
	{
		N_tk = mxGetN(prhs[0]);
	}

//	printf("N_tk = %d\n",N_tk);
	
	//plhs[0] = mxCreateDoubleMatrix(N_x,1, mxREAL);
    plhs[0] = mxCreateDoubleMatrix(mrows,ncols, mxREAL);
	
	y = mxGetPr(plhs[0]);
	tk = mxGetPr(prhs[0]);
	x = mxGetPr(prhs[1]);
    	
	for(ii=0;ii<N_x;ii++)
	{
        for(jj=0;jj<N_tk;jj++)
        {
  //          printf("x[ii] = %f, tk[jj] = %f\n",x[ii],tk[jj]);
            if( !( x[ii] > tk[jj]))
            {
        		y[ii] = (double)(jj);    
                break;  
            }
                
        }//for jj
        if( jj == N_tk)
        {
            y[ii] = N_tk;
        }
        
//        printf("y[ii] = %f\n",y[ii]);
//        printf("---------\n");

	}//for ii
	
    return;
}
