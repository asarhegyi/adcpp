#include "math.h"
#include "mex.h"   //--This one is required

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double  *x;//input
	double *y;//output
	
     //mwSize mrows,ncols;
    unsigned int mrows,ncols;
	double *max_value;//tranzition levels
	
	unsigned int N_x,M_max;
	unsigned int temp_index;
	
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
		mexErrMsgTxt("Data must be a vector.");
	}

	if( ( (mxGetN(prhs[1]) != 1) || (mxGetM(prhs[1]) != 1)))
	{
		mexErrMsgTxt("Maximum value must be a scalar.");
	}
	
	if(mxGetN(prhs[0]) == 1) 
	{
		N_x = mxGetM(prhs[0]);
	}
	else
	{
		N_x = mxGetN(prhs[0]);
	}

	max_value = mxGetPr(prhs[1]);
	M_max = *max_value+1;
	
    plhs[0] = mxCreateDoubleMatrix(M_max,1, mxREAL);
	
	y = mxGetPr(plhs[0]);
	x = mxGetPr(prhs[0]);
	
	for(jj=0;jj<M_max;jj++)
	{
			y[jj] = 0;
	}//for jj
	
	for(ii=0;ii<N_x;ii++)
	{
		if( (x[ii] <= M_max) && (0.0<=x[ii])) 
		{
			temp_index = (unsigned int)x[ii];
			//printf("temp_index = %d\n",temp_index);
			y[temp_index] += 1;
		}//if

	}//for ii
	
    return;
}
