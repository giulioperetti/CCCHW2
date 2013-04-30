// BINARY CONVOLUTIONAL ENCODER

#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // input output check
    if (nrhs != 3)
        mexErrMsgTxt("Three input arguments required");
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required");
    
    // declare variables
    mxArray *u_in_m, *S_in_m, *O_in_m, *y_out_m;
    const mwSize *u_dims, *S_dims, *O_dims, *y_dims;
    double *u, *S, *O, *y;
    int mu, numdims;
    int i;
    
    // associate inputs
    u_in_m = mxDuplicateArray(prhs[0]);
    S_in_m = mxDuplicateArray(prhs[1]);
    O_in_m = mxDuplicateArray(prhs[2]);
    
    // figure out dimensions
    for (i=0;i<3;i++)
    {
        numdims = mxGetNumberOfDimensions(prhs[i]);
        if (numdims != 2)
            mexErrMsgTxt("Only 2 dimensions input matrices allowed");
    }
    
    u_dims = mxGetDimensions(prhs[0]);
    S_dims = mxGetDimensions(prhs[1]);
    O_dims = mxGetDimensions(prhs[2]);
   
    
    mu = (int)u_dims[1]; // input length
        
    // associate outputs
    y_out_m = plhs[0] = mxCreateDoubleMatrix(u_dims[0],u_dims[1],mxREAL);
    
    // associate pointers
    u = mxGetPr(u_in_m);
    S = mxGetPr(S_in_m);
    O = mxGetPr(O_in_m);
    y = mxGetPr(y_out_m);
    
    double s = 0;
    
    for(i=0;i<mu;i++)
    {       
        y[i] = O[(int)s+(int)u[i]*O_dims[0]];    
        s = S[(int)s+(int)u[i]*S_dims[0]];
       
    }
    
    return;
 }