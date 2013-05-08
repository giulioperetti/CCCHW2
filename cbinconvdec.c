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
    
void stamp(double **M,int r, int c);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // input output arguments check
    if (nrhs != 4)
        mexErrMsgTxt("Four input arguments required");
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required");
    
    // variables declaration
    mxArray *r_in_m, *S_in_m, *O_in_m, *N_in_m, *u_hat_out_m;
    const mwSize *r_dims, *S_dims, *O_dims, *N_dims, *u_hat_dims;
    double **r, *rmem, *u_hat, **S, *Smem, **O, *Omem, **N, *Nmem, *m, *C, *C_aux, **U, *Umem, **U_aux, *U_auxmem;
    int *p, t, max, mu, nu, n_states, i, j, k;
    
    // inputs association
    r_in_m = mxDuplicateArray(prhs[0]);
    S_in_m = mxDuplicateArray(prhs[1]);
    O_in_m = mxDuplicateArray(prhs[2]);
    N_in_m = mxDuplicateArray(prhs[3]);
    
    // dimensions check
    for (i=0;i<4;i++)
        if (mxGetNumberOfDimensions(prhs[i]) != 2)
            mexErrMsgTxt("Only 2 dimensions input matrices allowed");
    
    r_dims = mxGetDimensions(prhs[0]);  
    S_dims = mxGetDimensions(prhs[1]);
    O_dims = mxGetDimensions(prhs[2]);
    N_dims = mxGetDimensions(prhs[3]);
    if (r_dims[1] != 2)
        mexErrMsgTxt("The first input argument must have 2 columns");
    if (S_dims[0] != 2)
        mexErrMsgTxt("The fourth input argument must have 2 rows");
    if (O_dims[0] != 2)
        mexErrMsgTxt("The fourth input argument must have 2 rows");
    if (N_dims[0] != 2)
        mexErrMsgTxt("The fourth input argument must have 2 rows");
    
    // input length, number of states, system memory
    mu = (int)r_dims[0]; 
    n_states = S_dims[1]; 
    nu = (int)ceil(log(n_states)/log(2)); 
    
    // outputs association
    u_hat_out_m = plhs[0] = mxCreateDoubleMatrix(1,mu,mxREAL);
    
    // pointers association and allocation
    rmem = mxGetPr(r_in_m);
    r = mxCalloc(2, sizeof(double));
    r[0] = &rmem[0];
    r[1] = &rmem[mu];
        
    u_hat = mxGetPr(u_hat_out_m);
    
    Smem = mxGetPr(S_in_m);
    S = mxCalloc(n_states, sizeof(double));
    Omem = mxGetPr(O_in_m);
    O = mxCalloc(n_states, sizeof(double));
    Nmem = mxGetPr(N_in_m);
    N = mxCalloc(n_states, sizeof(double));
    for (i=0;i<n_states;i++)
    {
        S[i] = &Smem[i*2];
        O[i] = &Omem[i*2];
        N[i] = &Nmem[i*2];       
    } 
    
    // VITERBI DECODER
      
    // cost vectors association and allocation
    C = mxCalloc(n_states, sizeof(double)); 
    C_aux = mxCalloc(n_states, sizeof(double)); 
    
    for (i=0;i<n_states;i++)
    {
        C[i] = -DBL_MAX;
        C_aux[i] = 0;
    }
    
    C[0] = 0;
     
    // survivors matrix association and allocation
    Umem = mxCalloc(n_states*mu, sizeof(double));
    U = mxCalloc(n_states, sizeof(double));
    U_auxmem = mxCalloc(n_states*(15*nu), sizeof(double));
    U_aux = mxCalloc(n_states, sizeof(double));
    
    for (i=0;i<n_states;i++)
    {
        U[i] = &Umem[i*mu];
        U_aux[i] = &U_auxmem[i*15*nu];
    }
    
    // auxiliary variables association and allocation
    p = mxCalloc(2,sizeof(int));
    m = mxCalloc(2,sizeof(double));
    
    
    // Viterbi algorithm
    for (i=0;i<mu;i++)
    {
        for (j=0;j<n_states;j++)
        {  
            // predecessors
            p[0] = (int)N[j][0];
            p[1] = (int)N[j][1];
                        
            // transition   
            t = floor(j/2);
            
            // paths metric
            m[0] = C[p[0]] + r[0][i]*pamap(floor(O[p[0]][t]/2))
                           + r[1][i]*pamap(floor((int)O[p[0]][t]%2));
            
            m[1] = C[p[1]] + r[0][i]*pamap(floor(O[p[1]][t]/2))
                           + r[1][i]*pamap(floor((int)O[p[1]][t]%2));
                        
            // maximum metric selection
            if (m[0] > m[1])
                max = 0;
            else max = 1;

            // (auxiliary) path and cost update
            C_aux[j] = m[max];

            for (k=0;k<min(i,15*nu);k++)
                U_aux[j][k] = U[p[max]][i-min(i,15*nu-1)+k];
            
            U_aux[j][min(i,15*nu-1)] = t;
        }
        
        // path and cost update
        max = C_aux[0];
        for (k=1;k<n_states;k++)
        {
            if (C_aux[k]<max) 
                    max = C_aux[k]; 
        }
        for (k=0;k<n_states;k++)
            C[k] = C_aux[k] - max;
        
        for (k=0;k<min(i+1,15*nu);k++)
            for (j=0;j<n_states;j++)
                U[j][max(0,i-15*nu+1)+k] = U_aux[j][k];        
         
    }
    
    // winner path selection
    max = 0;
    for (k=1;k<n_states;k++)
    {
        if (C[k]>C[max]) 
                max = k; 
    }
    
    // output 
    for (k=0;k<mu;k++)
    {
        u_hat[k] = U[max][k];
    }
    
    // disallocation
    mxFree(m);
    mxFree(p);
    mxFree(C);
    mxFree(C_aux);
    mxFree(Umem);
    mxFree(U);
    mxFree(U_aux);
    mxFree(U_auxmem);
    
    return;
}

int pamap(int x)
{
    // pam map: 0 -> -1 ; 1 -> 1
    return 2*x-1;
}