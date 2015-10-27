#include "Gauss.h"

bool SLAE_solution_Gauss(double **A, double *b, double *x, int n){
	bool bad;
	for(int i = 0 ; i < n; i++)
		x[i] = b[i];
	bad = trianglematrix1(A,x,n);
	if(!bad) return bad;
	for (int i = n-1; i >=0 ; i--){
		double sum = 0;
		for(int j = i+1; j < n; j++)
			sum += A[i][j]*x[j];
		x[i] -= sum;
	}
	return true;
}


bool trianglematrix1(double **A, double *x, int n){
	for(int i = 0 ; i < n ; i++){
		transf1(A,x,i,n);
		double A_d = A[i][i];
		if(fabs(A_d) < 1E-20) return false;

		for(int p = i; p < n; p++){
			A[i][p] /= A_d;
		}
		x[i] /= A_d;

		for(int j = i + 1; j < n; j++){
			double A_j = A[j][i];
			if(A_j != 0){
				for(int k = i; k < n; k++)
					A[j][k] -= A[i][k]*A_j;
				x[j] -= x[i]*A_j;

			}
		}
	}
	return true;
}


void transf1(double **A, double *x, int i, int n){

	int line = i;
	for(int j = i+1; j < n; j++)
		if(fabs(A[j][i]) > fabs(A[line][i]))
			line = j;

	if (line != i){
		double c;
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


