#pragma once
#include <math.h>
#include <stdio.h>
class MSG_my
{
	public:
		void give_data(int *ig1,int *jg1,double *gg1, double *di1, double *pr1, double *x01, int m1, int ngg1);	

	private:
		int *ig, *jg, m, n_gg;
		int time;
		double *gg, *di, *pr, *x0, *q2;
		double *Ll;
		double *Ld;
		bool flag;
		void MSG_X();
		void MSG_DI();
		double skal(double *a, double *b);
		void mult_matrix(double *x, double *&b);
		void slau_X(double *&a, double *b);
		void slau_DI(double *&a, double *b);
		void LLT();
		void DI();
};