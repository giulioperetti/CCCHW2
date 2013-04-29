#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //expecting 
    if (nrhs != 1)
        mexErrMsgTxt("One input argument required");
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required");
    
    //declare variables
    mxArray *a_in_m, *b_out_m;
    const mwSize *dims;
    double *a, *b;
    int dimx, dimy, numdims;
    int i,j;
    
    //associate inputs
    a_in_m = mxDuplicateArray(prhs[0]);
    
    //figure out dimensions
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimy = (int)dims[0]; dimx = (int)dims[1];
    
    //associate outputs
    b_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
    
    //associate pointers
    a = mxGetPr(a_in_m);
    b = mxGetPr(b_out_m);
    
    //do something
    for(i=0;i<dimx;i++)
    {
        for(j=0;j<dimy;j++)
        {
            mexPrintf("element[%d][%d] = %f\n",j,i,a[i*dimy+j]);
            b[i*dimy+j] = a[i*dimy+j];
        }
    }
    
    return;
}