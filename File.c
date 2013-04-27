#include "math.h"
#include "complex.h"
#include "stdlib.h"
#include "mex.h"

//for i=1:mu
//    z = [u(i)+z(1:1+nud-1)*d(2:end).' z(1:end-1)];
//    y(i) = z(1:num+1)*n(1:num+1).';
//end
//y = mexfunction(u,z,d,num)

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    double *xR, *xI, *hR, *hI, *yR, *yI;
    
    int i, Nx, Nh;
    
    if (nrhs != 2)
        mexErrMsgTxt("Two input arguments required");
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required");
    
    Nx = mxGetNumberOfElements(prhs[0]);
    xR = mxGetPr(prhs[0]);
    xI = mxGetPi(prhs[0]);
    if (!xI)
        xI = mxCalloc(Nx, sizeof(double));
    
    Nh = mxGetNumberOfElements(prhs[1]);
    hR = mxGetPr(prhs[1]);
    hI = mxGetPi(prhs[1]);
    if (!hI)
        hI = mxCalloc(Nx, sizeof(double));
    
    plhs[0] = mxCreateDoubleMatrix(1, Nx, mxCOMPLEX);
    yR = mxGetPr(plhs[0]);
    yI = mxGetPi(plhs[0]);
    
    return;
    
}