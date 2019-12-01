#include "../include/for_you_to_do.h"
/**
 * 
 * this function computes LU factorization
 * for a square matrix
 * 
 * syntax 
 *  
 *  input : 
 *      A     n by n , square matrix
 *      ipiv  1 by n , vector
 *      n            , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf(double *A, int *ipiv, int n) 
{
    /* add your code here */
	int i,i1,j,k,t,maxind,temps; 
	double max_1,tempv;
    /* Factorize A */
	for(i=0; i<n-1; i++){
		maxind = i;
		max_1 = fabs(A[i*n+i]);
		for(t=i+1; t<n; t++){
			if(fabs(A[t*n+i])>max_1){
				maxind = t;
				max_1 = fabs(A[t*n+i]);
			}
		}
		if(max_1 == 0){
			return -1;
		}
		else{
			if(maxind != i){
				/* save pivoting information */
				temps = ipiv[i];
				ipiv[i] = ipiv[maxind];
				ipiv[maxind] = temps;
				
				/* swap rows */
				for(i1 = 0; i1 < n; i1++ ){
					tempv = A[i*n+i1];
					A[i*n+i1] = A[maxind*n+i1];
					A[maxind*n+i1] = tempv;
				}	
			}
		}
		/* factorization */
		for(j = i+1; j<n; j++){
			A[j*n+i] = A[j*n+i]/A[i*n+i];
			for(k=i+1; k<n; k++){
				A[j*n+k] = A[j*n+k] - A[j*n+i]*A[i*n+k];
			}
		}
	}
    return 0;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix . according to lecture slides, this
 * function computes forward AND backward subtitution in the
 * same function.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      none
 * 
 **/
void mydtrsv(char UPLO, double *A, double *B, int n, int *ipiv)
{
    /* add your code here */
    double y[n];
    /*double b[n];*/
    double x[n];
    int i,j;
    double r;
	double s;
    
    /* forward substitution */
    if(UPLO=='U'){
    	y[0] = B[ipiv[0]];
	    for(i = 1; i < n; i++){
	    	r = 0;
	    	for(j=0; j<i-1; j++){
				r += y[j]*A[i*n+j];
			}
	    	y[i] = B[ipiv[i]] - r;
		}
	}
	
	/* back substitution */
	if(UPLO=='L'){
		x[n-1] = y[n-1]/A[(n-1)*n+n-1];
		for(i = n-2; i>=0; i--){
			s = 0;
			for(j=i+1 ; j<n; j++){
				s += x[j]*A[i*n+j];
			}
			x[i] = (y[i]-s)/A[i*n+i];
		}
	}
	
    
    return;
}

/**
 * 
 * Same function as what you used in lab1, cache_part4.c : optimal( ... ).
 * 
 **/
void mydgemm(double *A, double *B, double *C, int n, int i, int j, int k, int b)
{
    /* add your code here */
    /* please just copy from your lab1 function optimal( ... ) */
	int i1,j1,k1;
	for (k = 0; k < n; k+=b)
		for (i = 0; i < n; i+=b)
			for (j = 0; j < n; j+=b)
				/* B x B mini matrix multiplications */
				for (i1 = i; i1 < i+b&&i1<n; i1+=3){
					
					for (j1 = j; j1 < j+b&&j1<n; j1+=3) {
							register int t = i1*n+j1;
							register int tt = t+n;
							register int ttt = tt+n;
	
							register double c00 = C[t];
							register double c01 = C[t+1];
							register double c02 = C[t+2];
							register double c10 = C[tt];
							register double c11 = C[tt+1];
							register double c12 = C[tt+2];
							register double c20 = C[ttt];
							register double c21 = C[ttt+1];
							register double c22 = C[ttt+2];
							
							
						for (k1 = k; k1 < k+b&&k1<n; k1+=3){
							register int ta = i1*n+k1;
							register int tta = ta+n;
							register int ttta = tta+n;
							
							register double a00 = A[ta];
							register double a10 = A[tta];
							register double a20 = A[ttta];
							
							register int tb = k1*n+j1;
							register int ttb = tb+n;
							register int tttb = ttb+n;

							register double b00 = B[tb];
							register double b01 = B[tb+1];
							register double b02 = B[tb+2];
							
							c00 += a00*b00;
							c01 += a00*b01;
							c02 += a00*b02;
							c10 += a10*b00;
							c11 += a10*b01;
							c12 += a10*b02;
							c20 += a20*b00;
							c21 += a20*b01;
							c22 += a20*b02;
							
							a00 = A[ta+1];
							a10 = A[tta+1];
							a20 = A[ttta+1];
							
							b00 = B[ttb];
							b01 = B[ttb+1];
							b02 = B[ttb+2];
							
							c00 += a00*b00;
							c01 += a00*b01;
							c02 += a00*b02;
							c10 += a10*b00;
							c11 += a10*b01;
							c12 += a10*b02;
							c20 += a20*b00;
							c21 += a20*b01;
							c22 += a20*b02; 
							
							a00 = A[ta+2];
							a10 = A[tta+2];
							a20 = A[ttta+2];
							
							b00 = B[tttb];
							b01 = B[tttb+1];
							b02 = B[tttb+2];
							
							c00 += a00*b00;
							c01 += a00*b01;
							c02 += a00*b02;
							c10 += a10*b00;
							c11 += a10*b01;
							c12 += a10*b02;
							c20 += a20*b00;
							c21 += a20*b01;
							c22 += a20*b02; 
	
					}	
							C[t] = c00;
							C[t+1] = c01;
							C[t+2] = c02;
							C[tt] = c10;
							C[tt+1] = c11;
							C[tt+2] = c12;
							C[ttt] = c20;
							C[ttt+1] = c21;
							C[ttt+2] =c22;	
						}
				}
    return;
}

/**
 * 
 * this function computes triangular matrix - vector solver
 * for a square matrix using block gepp introduced in course
 * lecture .
 * 
 * just implement the block algorithm you learned in class.
 * 
 * syntax 
 *  
 *  input :
 *      UPLO  'L' or 'U' , denotes whether input matrix is upper
 *                         lower triangular . ( forward / backward
 *                         substitution )
 * 
 *      A     n by n     , square matrix
 * 
 *      B     1 by n     , vector
 * 
 *      ipiv  1 by n     , vector , denotes interchanged index due
 *                                  to pivoting by mydgetrf()
 * 
 *      n                , length of vector / size of matrix
 *  
 *  output :
 *      return -1 : if the matrix A is singular (max pivot == 0)
 *      return  0 : return normally 
 * 
 **/
int mydgetrf_block(double *A, int *ipiv, int n, int b) 
{
	int i,i1,j,k,ib,max_1,maxind,temps;
	int row_start, col_start;	
	double max_1,tempv;
	
	for(row_start = 0; row_start < n; row_start+=b){
		col_start = row_start;
			if(n - row_start > 0){
				/*Factorize left blue*/ 
				for(i=row_start; i<n-1; i++){
					maxind = i;
					max_1 = fabs(A[i*n+i]);
					for(t=i+1; t<n; t++){
						if(fabs(A[t*n+i])>max_1){
							maxind = t;
							max_1 = fabs(A[t*n+i]);
						}
					}
					if(max_1 == 0){
						return -1;
					}
					else{
						if(maxind != i){
							/* save pivoting information */
							temps = ipiv[i];
							ipiv[i] = ipiv[maxind];
							ipiv[maxind] = temps;
							
							/* swap rows */
							for(i1 = col_start; i1 < n; i1++ ){
								tempv = A[i*n+i1];
								A[i*n+i1] = A[maxind*n+i1];
								A[maxind*n+i1] = tempv;
							}	
						}
					}
					/* factorization */
					for(j = i+1; j<n; j++){
						A[j*n+i] = A[j*n+i]/A[i*n+i];
						for(k=i+1; k<n; k++){
							A[j*n+k] = A[j*n+k] - A[j*n+i]*A[i*n+k];
						}
					}
				}
				/*End factorize left blue*/
				
				
				/*Calculate upper pink*/
				int l_row = row_start;
				int l_col = col_start;
				int rup_row = row_start;
				int rup_col = col_start;
				int ii;
				double sum_al=0;
				

				for(ii = rup_col+b; ii <  n; ii++){
					A[rup_row*n+ii] = A[rup_row*n+ii]/1;
				}
				
				for(rup_row = row_start+1; rup_row < rup_row + b&&row_start+b < n; rup_row++){
					for(rup_col = col_start+b; rup_col < n; rup_col++){
						sum_al = 0;
						for(ii = col_start; ii < rup_row; ii++){
							sum_al += A[rup_row*n+ii]*A[rup_col*n+ii];
						}
						A[rup_row*n+rup_col] = (A[rup_row*n+rup_col] - sum_al)/1;
					}
				}
				
				/*End calculate upper pink*/
				

				/*Green update*/
	
				  int j1,k1;
 			
                                /* B x B mini matrix multiplications */

                		     for (i = row_start+b; i < n; i+=b){
                        		    for (j = col_start+b; j < n; j+=b){
                                /* B x B mini matrix multiplications */
                                		for (i1 = i; i1 < i+b&&i1<n; i1+=3){

                                       			for (j1 = j; j1 < j+b&&j1<n; j1+=3) {       		

							register int t = i1*n+j1;
							register int tt = t+n;
							register int ttt = tt+n;
	
							register double c00 = A[t];
							register double c01 = A[t+1];
							register double c02 = A[t+2];
							register double c10 = A[tt];
							register double c11 = A[tt+1];
							register double c12 = A[tt+2];
							register double c20 = A[ttt];
							register double c21 = A[ttt+1];
							register double c22 = A[ttt+2];
							
						for (k1 = col_start; k1 < col_start+b ; k1+=3){
							register int ta = i1*n+k1;
							register int tta = ta+n;
							register int ttta = tta+n;
							
							register double a00 = A[ta];
							register double a10 = A[tta];
							register double a20 = A[ttta];
							
							register int tb = k1*n+j1;
							register int ttb = tb+n;
							register int tttb = ttb+n;

							register double b00 = A[tb];
							register double b01 = A[tb+1];
							register double b02 = A[tb+2];
							
							c00 -= a00*b00;
							c01 -= a00*b01;
							c02 -= a00*b02;
							c10 -= a10*b00;
							c11 -= a10*b01;
							c12 -= a10*b02;
							c20 -= a20*b00;
							c21 -= a20*b01;
							c22 -= a20*b02;
							
							a00 = A[ta+1];
							a10 = A[tta+1];
							a20 = A[ttta+1];
							
							b00 = A[ttb];
							b01 = A[ttb+1];
							b02 = A[ttb+2];
							
							c00 -= a00*b00;
							c01 -= a00*b01;
							c02 -= a00*b02;
							c10 -= a10*b00;
							c11 -= a10*b01;
							c12 -= a10*b02;
							c20 -= a20*b00;
							c21 -= a20*b01;
							c22 -= a20*b02; 
							
							a00 = A[ta+2];
							a10 = A[tta+2];
							a20 = A[ttta+2];
							
							b00 = A[tttb];
							b01 = A[tttb+1];
							b02 = A[tttb+2];
							
							c00 -= a00*b00;
							c01 -= a00*b01;
							c02 -= a00*b02;
							c10 -= a10*b00;
							c11 -= a10*b01;
							c12 -= a10*b02;
							c20 -= a20*b00;
							c21 -= a20*b01;
							c22 -= a20*b02; 

					}	
							A[t] = c00;
							A[t+1] = c01;
							A[t+2] = c02;
							A[tt] = c10;
							A[tt+1] = c11;
							A[tt+2] = c12;
							A[ttt] = c20;
							A[ttt+1] = c21;
							A[ttt+2] =c22;	
						}
				}	
				}
				}
				
				
				/*Green update*/
			}
		
		
	}
	
	

	
	
	
	
    return 0;
}



