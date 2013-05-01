// BINARY CONVOLUTIONAL DECODER

#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // input output arguments check
    if (nrhs != 4)
        mexErrMsgTxt("Four input arguments required");
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required");
    
    // declare variables
    mxArray *r_in_m, *S_in_m, *O_in_m, *N_in_m, *u_hat_out_m;
    const mwSize *r_dims, *S_dims, *O_dims, *N_dims, *u_hat_dims;
    double *r, *S, *O, *N, *u_hat, *C, *C_aux, *U, *U_aux, s;
    int *p, *t, mu, nu, i, j;
    
    // associate inputs
    r_in_m = mxDuplicateArray(prhs[0]);
    S_in_m = mxDuplicateArray(prhs[1]);
    O_in_m = mxDuplicateArray(prhs[2]);
    N_in_m = mxDuplicateArray(prhs[3]);
    
    // figure out dimensions
    for (i=0;i<4;i++)
    {
        if (mxGetNumberOfDimensions(prhs[i]) != 2)
            mexErrMsgTxt("Only 2 dimensions input matrices allowed");
    }
    
    r_dims = mxGetDimensions(prhs[0]);
    if (r_dims[0] != 2)
        mexErrMsgTxt("The first input argument must have 2 rows");
    S_dims = mxGetDimensions(prhs[1]);
    if (S_dims[1] != 2)
        mexErrMsgTxt("The second input argument must have 2 columns");
    O_dims = mxGetDimensions(prhs[2]);
    if (O_dims[1] != 2)
        mexErrMsgTxt("The third input argument must have 2 columns");
    N_dims = mxGetDimensions(prhs[3]);
    if (N_dims[1] != 2)
        mexErrMsgTxt("The fourth input argument must have 2 columns");
    
    mu = (int)r_dims[1]; // input length
    nu = (int)ceil(log(S_dims[0])/log(2)); // system memory 
    
    // associate outputs
    u_hat_out_m = plhs[0] = mxCreateDoubleMatrix(1,mu,mxREAL);
    
    // associate pointers
    r = mxGetPr(r_in_m);
    S = mxGetPr(S_in_m);
    O = mxGetPr(O_in_m);
    N = mxGetPr(N_in_m);
    u_hat = mxGetPr(u_hat_out_m);
    
    /* VITERBI DECODING */
      
    // cost vectors initialization
    C = mxCalloc(S_dims[0], sizeof(double)); 
    C_aux = mxCalloc(S_dims[0], sizeof(double)); 
    
    for (i=0;i<S_dims[0];i++)
    {
        C[i] = -1000;
        C_aux[i] = 0;
    }
    
    C[0] = 0;
     
    // survivors matrix initialization
    
    U = mxCalloc(S_dims[0]*mu, sizeof(double));
    U_aux = mxCalloc(S_dims[0]*(5*nu), sizeof(double));
    
    // state and auxiliary variables initialization
    s = 0;
    p = mxCalloc(2,sizeof(int));
    t = mxCalloc(2,sizeof(int));
    
    // Viterbi algorithm
    
    for (i=0;i<mu;i++)
    {
        for (j=0;j<S_dims[0];j++)
        {
            // analyzing state j
           
            // predecessors
            p[0] = (int)N[j];
            p[1] = (int)N[j+S_dims[0]];
            
            // transitions
            
            if (S[p[0]] == j)
                t[0] = 0;
            else t[0] = 1;
            if (S[p[1]] == j)
                t[1] = 0;
            else t[1] = 1;
            
            mexPrintf("State %d: p = [%d][%d]",j,p[0],p[1]);
//             % transitions from predecessors
//             trans = [find(S(pred(1)+1,:)==j-1)-1 
//                 find(S(pred(2)+1,:)==j-1)-1];
//             
//             % cost function
//             temp = [C(pred(1)+1) + r(:,i).'*pamap(de2bi(O(pred(1)+1,trans(1)+1),2,'left-msb')).' 
//                 C(pred(2)+1) + r(:,i).'*pamap(de2bi(O(pred(2)+1,trans(2)+1),2,'left-msb')).'];
//             
//             % select the max value
//             [mcost ind] = max(temp);
//             
//             % path and cost update
//             Caux(j) = mcost;
//             Uaux(j,1:min(i,5*nu)) = U(pred(ind)+1,max(i-5*nu+1,1):i);
//             Uaux(j,min(i,5*nu)) = trans(ind);
        }
    }
    
//     for(i=0;i<mu;i++)
//     {       
//         yd = O[(int)s+(int)u[i]*O_dims[0]]; // yd is in {0,1,2,3}  
//         s = S[(int)s+(int)u[i]*S_dims[0]]; // state updating
//         y[2*i] = floor(yd / (double)2);
//         y[2*i+1] = (double)((int)yd % 2);
//     }
    
    return;
 }