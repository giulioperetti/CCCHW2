// binconvmaps

#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // input output check
    if (nrhs != 1)
        mexErrMsgTxt("One input argument required");
//    if (nlhs != 3)
//        mexErrMsgTxt("Three output arguments required");
    
    // declare variables
    mxArray *n_in_m, *S_out_m, *O_out_m, *N_out_m;
    const mwSize *dims;
    double *n, *S;
    int dimx, dimy, numdims;
    int i,j;
    
    // associate inputs
    n_in_m = mxDuplicateArray(prhs[0]);
    
    // figure out dimensions
    numdims = mxGetNumberOfDimensions(prhs[0]);
    if (numdims != 2)
        mexErrMsgTxt("Only 2 dimensions input matrices allowed");
    dims = mxGetDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];
    
    
    
    // associate outputs
    S_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    
    // associate pointers
    n = mxGetPr(n_in_m);
    S = mxGetPr(S_out_m);
    
    for(i=0;i<dimx;i++)
    {
        for(j=0;j<dimy;j++)
        {
            mexPrintf("element[%d][%d] = %f\n",j,i,n[i*dimy+j]);
            S[i*dimy+j] = n[i*dimy+j];
        }
    }
    
    //mexPrintf("numdims = %d",numdims);
    
    return;
 }