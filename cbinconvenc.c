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
    double *u, *S, *O, *y, yd;
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
    if (u_dims[0] != 1)
        mexErrMsgTxt("The first input argument must have only one row");
    S_dims = mxGetDimensions(prhs[1]);
    if (S_dims[1] != 2)
        mexErrMsgTxt("The second input argument must have 2 columns");
    O_dims = mxGetDimensions(prhs[2]);
    if (O_dims[1] != 2)
        mexErrMsgTxt("The third input argument must have 2 columns");
    
    mu = (int)u_dims[1];
        
    // associate outputs
    y_out_m = plhs[0] = mxCreateDoubleMatrix(2,mu,mxREAL);
    
    // associate pointers
    u = mxGetPr(u_in_m);
    S = mxGetPr(S_in_m);
    O = mxGetPr(O_in_m);
    y = mxGetPr(y_out_m);
    
    
    // code generation
    double s = 0;
    
    for(i=0;i<mu;i++)
    {       
        // decimal representation of y[i]
        yd = O[(int)s+(int)u[i]*O_dims[0]];
        
        // state updating
        s = S[(int)s+(int)u[i]*S_dims[0]]; 
        
        // decimal to binary conversion
        y[2*i] = floor(yd / (double)2);
        y[2*i+1] = (double)((int)yd % 2);
    }
    
    return;
 }