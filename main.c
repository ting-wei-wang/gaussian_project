/*************************************************************
* In this sample program. the power method and inverse power
* method are used to compute the max and min eigenvalues of
* a matrix.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "definition.h"
#include "vec_mtx.h"

/*--- Declare the coef. matrix, the unkown vec. and the rhs.. ---*/
double   **A, *x, *b , **At ,**all,*vec_new,big,total,y,v;
int      n,m,i; // dimension of the system.




/*------------------------------------------------------
 * Create a Hilbert linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix.
 */
void tranpose_matrix (double **A, double *b, int n,int m){
	int i ,j;
	double num = 1.0;
	for(i = 0 ;i < n ;i++){
		for(j = 0; j < m ;j++){
			A[j][i] = pow( num ,j );
			
		}
		
		num = num + 0.2;
	}
}
void Hilbert_linear_system(double **A, double *b, int n,int m)
{
	int  i, j;
	double num = 1.0;
	double sum ;
	double no;
    
	for(i = 0;i< n ; i++){
			A[i][0] = 1;
			//printf("%lf ",A[i][0]);
		}
	for(i = 0;i< n ; i++){
		//for(num = 1 ; num <= 3; num = num + 0.2)
		sum = 1;
			for(j = 1 ; j < m ;j++ ){
			// A[i][j] = pow(num ,j);
				
				sum = sum*num;
				
				A[i][j] = sum;
			}
			//printf("%lf ",A[i][1]);
			num = num + 0.2;
			
		}
	for(i = 0 ;i < n;i++){
		no = 0;
		for(j = 0;j < m ;j++){
			no = no + A[i][j];
		}
		b[i] = no;
		}
	
}




/*------------------------------------------------------
 * Create a symmetric linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix.
 */
void symmetric_linear_system(double **A, double *b, int n)
{
	int  i, j;

	for(i=0;i<n;i++){
		for(j=i;j<n;j++){
			A[j][i] = rand()%10 +1.0;
			A[i][j] = rand() % 10 + 1.0;
		}
	}
	for(i=0;i<n;i++){
		b[i] = 0.0;
		for(j=0;j<n;j++) b[i] += A[i][j];
	}
}

double Horner( double x, int m)
{
  double sum;
  int    i;

  sum = 1;
  for(i=0;i< m-1;i++){
    sum = sum*(x) +1;
  }
  return sum;
}

/*----------------------------------------------------
 * The main procedure
 */
int main(int argc, char **argv)
{

	n = 11;
	m = 8;
	
    A = alloc_mtx(n , m);
	At =  alloc_mtx(m , n);
    x = alloc_vec(m );
    b = alloc_vec(n);
	all = alloc_mtx(m,m);
	vec_new = alloc_vec(m);

    Hilbert_linear_system(A, b, n,m);
	fprintf(stderr,"A[][]=\n");
    print_mtx(A, n ,m);
	tranpose_matrix (At, b, n,m);
	//symmetric_linear_system(A, b, n);
// Print out the initial linear system
	
    fprintf(stderr,"A^T[][]=\n");
    print_mtx(At, m ,n);

	mtx_mtx_mult(all, At, A, m, n);
	fprintf(stderr,"ALL[][]=\n");
	 print_mtx(all, m ,m);

	 fprintf(stderr,"Horner:\n");
	 for(v = 1.0;v < 3.1 ; v = v + 0.2){

        y = Horner(v,m);
		printf(" %lf ",y);
	 }
	 printf("\n");
	fprintf(stderr,"b[]=\n");
	print_vec(b, n);

	 mtx_vec_mult(vec_new, At, b, m,n);
	fprintf(stderr,"A^T * b[]=\n");
	print_vec(vec_new, m);
// Perform Gaussian elimination
    gauss_elm(all, vec_new, m);
	//printf("--------------------------");
	//print_vec(vec_new,15 );
	fprintf(stderr," After forward elimination, A[][]=\n");
    print_mtx(all, m,m);
	print_vec(vec_new, m);
	//Perform backward substitution.
    back_substitute(all, x, vec_new, m);
    //print out the results
	fprintf(stderr,"The solution x[]=\n");
	print_vec(x, m);
	
	printf("infinite norm:\n");
	big = fabs(x[0] - 1);
	for(i = 1 ; i < m;i++){
      if(fabs(x[i] - 1) >= big){
		  big = fabs(x[i] - 1);
	}
	}
	printf("%lf\n",big);
	printf("----------------------\n");
	printf("two norm:\n");
	total = 0;
	for(i = 0 ; i < m;i++){
      total += ((x[i]-1)*(x[i]-1));
	}
	total = sqrt(total);
	printf("%lf\n",total);
	printf("----------------------------\n");
	
	system("pause");
	return 0;
}
