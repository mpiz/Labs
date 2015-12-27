#include "Gauss.h"

bool SLAE_solution_Gauss(dyn_matrix& A, dcomplex *b, dcomplex *x, int n){
	bool bad;
	for(int i = 0 ; i < n; i++)
		x[i] = b[i];
	bad = trianglematrix1(A,x,n);
	if(!bad) return bad;
	for (int i = n-1; i >=0 ; i--){
		dcomplex sum = 0.0;
		for(int j = i+1; j < n; j++)
			sum += A[i][j]*x[j];
		x[i] -= sum;
	}
	return true;
}


bool trianglematrix1(dyn_matrix& A, dcomplex *x, int n){
	for(int i = 0 ; i < n ; i++){
		transf1(A,x,i,n);
		dcomplex A_d = A[i][i];
		if(abs(A_d) < 1E-20) return false;

		for(int p = i; p < n; p++){
			A[i][p] /= A_d;
		}
		x[i] /= A_d;

		for(int j = i + 1; j < n; j++){
			dcomplex A_j = A[j][i];
			if(A_j != dcomplex(0, 0)){
				for(int k = i; k < n; k++)
					A[j][k] -= A[i][k]*A_j;
				x[j] -= x[i]*A_j;

			}
		}
	}
	return true;
}


void transf1(dyn_matrix& A, dcomplex *x, int i, int n){

	int line = i;
	for(int j = i+1; j < n; j++)
		if(abs(A[j][i]) > abs(A[line][i]))
			line = j;

	if (line != i){
		dcomplex c;
		for(int j = i; j < n; j++){
			c = A[i][j];
			A[i][j] = A[line][j];
			A[line][j] = c;
		}
		c = x[i];
		x[i] = x[line];
		x[line] = c;
	}
}


