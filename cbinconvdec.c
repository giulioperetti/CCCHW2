// BINARY CONVOLUTIONAL DECODER

#include <math.h>
#include <matrix.h>
#include <mex.h>
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif
    
int pamap(int x);
    
void stamp(double *M,int r, int c);

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
    int *p, *t, *m, max, mu, nu, i, j, k;
    
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
    m = mxCalloc(2,sizeof(int));
    
    // Viterbi algorithm
    for (i=0;i<mu;i++)
    {
        mexPrintf("\nSYMBOL %d:\n",i);
        for (j=0;j<S_dims[0];j++)
        {  
            mexPrintf(" State %d:\n",j);
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
            
            // paths metric
            m[0] = C[p[0]] + r[i*2]*pamap(floor(O[p[0]+t[0]*O_dims[0]]/2))
                           + r[i*2+1]*pamap((int)O[p[0]+t[0]*O_dims[0]]%2);
            
            mexPrintf(" C[p[0]] = %d",(int)C[p[0]]);
            mexPrintf(" r[i*2] = %d",(int)r[i*2]);
            mexPrintf(" r[i*2+1] = %d",(int)r[i*2+1]);
            mexPrintf(" PAM(0) = %d",(int)pamap(floor(O[p[0]+t[0]*O_dims[0]]/2)));
            
            m[1] = C[p[1]] + r[i*2]*pamap(floor(O[p[1]+t[1]*O_dims[0]]/2))
                           + r[i*2+1]*pamap((int)O[p[1]+t[1]*O_dims[0]]%2);
            mexPrintf("   Predecessors %d,%d; transitions %d,%d: metrics %d,%d\n",p[0],p[1],t[0],t[1],m[0],m[1]);
    
            // maximum metric selection
            if (m[0] > m[1])
                max = 0;
            else max = 1;
            
            mexPrintf("   Winner %d with metric %d and related transition %d\n",p[max],m[max],t[max]);
            
            // path and cost update
            C_aux[j] = m[max];
            mexPrintf("   C_aux[%d] = m[%d] = %d",j,max,m[max]);
            for (k=0;k<min(i,5*nu);k++)
            {
                mexPrintf("    U[%d] -> U_aux[%d] = %f\n",p[max]+S_dims[0]*(i-min(i,5*nu)+k),j+S_dims[0]*k,U_aux[j+S_dims[0]*k]);
                U_aux[j+S_dims[0]*k] = U[p[max]+S_dims[0]*(i-min(i,5*nu)+k)];
                
            }
            U_aux[j+S_dims[0]*min(i,5*nu-1)] = (double)t[max];
            mexPrintf("    t[%d] -> U_aux[%d] = %d\n",max,j+S_dims[0]*min(i,5*nu),t[max]);
            
            //stamp(U_aux,4,3);
            
        }
        
        //stamp(C_aux,4,1);
        
        // path and cost update
        max = C_aux[0];
        for (k=1;k<S_dims[0];k++)
        {
            if (C_aux[k]<max) 
                    max = C_aux[k]; 
        }
        for (k=0;k<S_dims[0];k++)
        {
            //C[k] = C_aux[k] - max;
            C[k] = C_aux[k];
        }
        for (k=0;k<S_dims[0]*min(i+1,5*nu);k++)
        {
            U[S_dims[0]*max(0,i-5*nu)+k] = U_aux[k];
           // mexPrintf("     U[%d] = U_aux[%d] = %f\n",S_dims[0]*max(0,i-5*nu)+k,k,U_aux[k]);
        }
//         C = Caux - max(Caux);
//         U(:,max(i-5*nu+1,1):i) = Uaux(:,1:min(i,5*nu));
         
    }
    
    max = 0;
    for (k=1;k<S_dims[0];k++)
    {
        if (C[k]>C[max]) 
                max = k; 
    }
    
    for (k=0;k<mu;k++)
    {
        u_hat[k] = U[max+k*S_dims[0]];
    }
    
    mxFree(m);
    mxFree(p);
    mxFree(t);
    mxFree(C);
    mxFree(C_aux);
    mxFree(U);
    mxFree(U_aux);
    
    return;
 }



void stamp(double *M,int r, int c)
{
    int i,j;
    for (i=0;i<r;i++)
    {
        for (j=0;j<c;j++)
        {
            mexPrintf("%d",(int)M[i+j*r]);
        }
        mexPrintf("\n");
    }
}

int pamap(int x)
{
    return 2*x-1;
}