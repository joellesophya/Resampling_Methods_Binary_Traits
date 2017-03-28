#include "svdcomp.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "nrutil.h"

void svdcomp(double **a, int m, int n, double w[], double **v) {
    /* a original matrix, a[1...m][1...n] */
    /* m number of rows of a */
    /* n number of columns of a */ 
    /* vector of length n */
    /* v is the V in a = U * M * V^T, v[1...n][1...n] */
    
    int i,j;
    
    /* create a dense matrix for the SVDLIBC and initalize with current a */
    double *A = dvector(0,m*n-1);
    
    for(j=1 ; j<=n ; j++) { /* ind j */
      for(i=1 ; i<=m ; i++) { /* ind i */
	A[(i-1) + (j-1)*m] = a[i][j];  
      }
    }
    
    long int M=m,N=n,LDA=m,LDU=m,LVT=n,LWORK=10*N,INFO,ONE=1;
    double *WORK = dvector(0,LWORK-1);
    double *S = dvector(0,N-1);
    double *VT = dvector(0,N*N-1);
    char cu = 'N', cv='A';
    dgesvd_(&cu, &cv, &M, &N, A, &LDA, S, 0, &ONE, VT, &N, WORK, &LWORK, &INFO);
 
    /* save the SVD results into the w vector and v matrix in correct form*/
    for(i=1 ; i<=n; i++) {
      w[i] = S[i-1];
    }
    
    for(i=1 ; i<=n ; i++) {
        for(j=1 ; j<=n ; j++) {
	  v[i][j] = VT[(j-1) + (i-1)*n];
        }
    }

    /* free the sparse matrix */
    free_dvector(A,0,m*n-1);
    free_dvector(WORK,0,LWORK-1);
    free_dvector(S,0,N-1);
    free_dvector(VT,0,N*N-1);
}



