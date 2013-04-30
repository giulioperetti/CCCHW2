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
    int nu, mu, dimx, dimy, numdims;
    int i,j;
    
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
    nu = ceil(log((int)S_dims[0])/log(2)); // system memory
    
    
   // if ((int)y_dims[1] != mu)
            //mexErrMsgTxt("Dimensions of u and y mismatch");
    
//     for (i=1;i<=mu;i++)
//         y(:,i) = de2bi(O(s+1,u(i)+1),2,'left-msb').';
//         s = S(s+1,u(i)+1);
//     end  
                
     dimy = (int)u_dims[0]; dimx = (int)u_dims[1];
//     
//     
//     
//     // associate outputs
//     S_out_m = plhs[0] = mxCreateDoubleMatrix(dimy,dimx,mxREAL);
//     
    // associate pointers
    u = mxGetPr(u_in_m);
    S = mxGetPr(S_in_m);
    O = mxGetPr(O_in_m);
    y = mxGetPr(y_out_m);
    
    int s = 0;
     
    for(i=0;i<dimx;i++)
    {       
            mexPrintf("y[%d] = %f\n",i,O[s*2+(int)u[i]]);
            
            s = S[s*2+i+1];
            mexPrintf("s = %d\n",s);
            //S[i*dimy+j] = n[i*dimy+j];
       
    }
    
   // mexPrintf("mu, (int)y_dims[1] = %d, %d",mu,(int)y_dims[1]);
    
    return;
 }