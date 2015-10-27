#include "solver_MSG.h"


void MSG_my::give_data(int *ig1,int *jg1,double *gg1, double *di1, double *pr1, double *x01, int m1, int n_gg1)
{
	ig = ig1;
	jg = jg1;
	gg = gg1;
	di = di1;
	pr = pr1;
	x0 = x01;
	
	m = m1;
	n_gg = n_gg1;
	
	MSG_X();
	
}

void MSG_my::MSG_X()
{

	double *z0 = new double[m];
	double *z1 = new double[m];
	double *r0 = new double[m];
	double *r1 = new double[m];
	double *xr = new double[m];
	double *xm = new double[m];
	double e = 1E-14;
	double maxiter = 1000000;
	double normaf=0;
	for(int j=0; j<m; j++)
		normaf+=pr[j]*pr[j];
	normaf = sqrt(normaf);

	Ll = new double[n_gg];
	Ld = new double[m];


	double a;
	double b1;
	double r;
	bool not_end = true;
	LLT();

	//невязка и подготовка к циклу
	mult_matrix(x0,xr);
	for(int i=0; i<m; i++)
	{
		r0[i] = pr[i] - xr[i];
	}
	slau_X(z0,r0);


	if (skal(r0,r0) == 0) not_end = false; 
	
	for(int i=1; i<maxiter && not_end; i++)
	{
		
		slau_X(xm,r0);	
		r = skal(xm,r0);

		mult_matrix(z0,xr);
		a = r/skal(xr,z0);

		for(int j=0; j<m; j++)
		{
			x0[j] = x0[j] + a*z0[j];
			r0[j] = r0[j] - a*xr[j];
		}

		slau_X(xm,r0);
		b1 = skal(xm,r0)/r;	


		for(int j=0; j<m; j++)
			z0[j] = xm[j] + b1*z0[j];

		r = skal(r0,r0);
		double normar = sqrt(r);	

		//if(i%100==0) printf("%e\n",normar/normaf);

		if((normar/normaf)<e) not_end = false;		

	}	
		

	delete[] z0;
	delete[] z1;
	delete[] r0;
	delete[] r1;
	delete[] xr;
	delete[] xm;
	delete[] Ll;
	delete[] Ld;

}
void MSG_my::MSG_DI()
{
	double *z0 = new double[m];
	double *z1 = new double[m];
	double *r0 = new double[m];
	double *r1 = new double[m];
	double *xr = new double[m];
	double *xm = new double[m];
	double e = 1E-16;
	double maxiter = 100000;
	double normaf=0;
	for(int j=0; j<m; j++)
		normaf+=pr[j]*pr[j];
	normaf = sqrt(normaf);

	Ll = new double[n_gg];
	Ld = new double[m];


	double a;
	double b1;
	double r;
	bool not_end = true;
	DI();

	//невязка и подготовка к циклу
	mult_matrix(x0,xr);
	for(int i=0; i<m; i++)
	{
		r0[i] = pr[i] - xr[i];
	}

	
	slau_DI(z0,r0);


	if (skal(r0,r0) == 0) not_end = false; 

	for(int i=1; i<maxiter && not_end; i++)
	{
		
		slau_DI(xm,r0);	
		r = skal(xm,r0);

		mult_matrix(z0,xr);
		a = r/skal(xr,z0);

		for(int j=0; j<m; j++)
		{
			x0[j] = x0[j] + a*z0[j];
			r0[j] = r0[j] - a*xr[j];
		}

		slau_DI(xm,r0);
		b1 = skal(xm,r0)/r;	


		for(int j=0; j<m; j++)
			z0[j] = xm[j] + b1*z0[j];

		r = skal(r0,r0);
		double normar = sqrt(r);	

		if(i%100==0) printf("%d %e \n", i, normar/normaf);

		if((normar/normaf)<e) not_end = false;		

	}	



	for(int i=0; i<m; i++)
		q2[i] = x0[i];

	delete[] z0;
	delete[] z1;
	delete[] r0;
	delete[] r1;
	delete[] xr;
	delete[] xm;
	delete[] Ll;
	delete[] Ld;

}
double MSG_my::skal(double *a, double *b)
{
	double s=0;
	for(int i=0; i<m; i++)
		s +=a[i]*b[i];
	return s;
}

void MSG_my::mult_matrix(double *x, double *&y)
{
	for(int i=0; i<m; i++)
	{
		y[i] = di[i]*x[i];		

		for(int j=ig[i]; j<ig[i+1]; j++)
		{			
			y[jg[j]] += gg[j]*x[i];
			y[i] += gg[j]*x[jg[j]];
		}
	}
}
void MSG_my :: slau_DI(double *&a, double *y)
{
	//прямой ход
	double *x = new double[m];
	for(int i=0; i<m; i++)
	{
		x[i] = (y[i])/Ld[i];
	}

	//обратный ход
		for(int k = m, k1 = m-1; k > 0; k--, k1--)
		{
			a[k1] = x[k1]/Ld[k1];				
		}
		delete[] x;

}

void MSG_my::slau_X(double *&a, double *y)
{
	//прямой ход
	double *x = new double[m];
	for(int i=0; i<m; i++)
	{
		int i0 = ig[i];
		int i1 = ig[i+1];
		
		double s = 0;
		for(int k=i0; k < i1; k++)
		{				
			s += Ll[k]*x[jg[k]];
		}

		x[i] = (y[i] - s)/Ld[i];
	}

	//обратный ход
		for(int k = m, k1 = m-1; k > 0; k--, k1--)
		{
			a[k1] = x[k1]/Ld[k1];
				for(int i = ig[k1]; i < ig[k]; i++)
					x[jg[i]] -= Ll[i]*a[k1];
		}
		delete[] x;

}
void MSG_my::LLT()
{
	//подготовка
	double sum_l = 0, sum_d = 0;
	for(int i=0; i<n_gg; i++)
		Ll[i] = gg[i];
			
	
		
	for(int i=0; i<m; i++)
		Ld[i] = di[i];
			

	//разложение
	for(int i=0; i<m; i++)
	{
		sum_d = 0;
		int i0 = ig[i];
		int i1 = ig[i+1];

		for(int k = i0; k < i1; k++)
		{
			sum_l = 0; 
			int j0 = ig[jg[k]];
			int j1 = ig[jg[k]+1];

			for(int l = i0; l < k; l++)
			{
				for(int n = j0; n < j1; n++)
				{
					if(jg[n] == jg[l])
					{
						sum_l += Ll[l]*Ll[n];
						j0++;
					}
				}
			}
			Ll[k] = (Ll[k] -  sum_l)/Ld[jg[k]];

			sum_d += Ll[k]*Ll[k];
		}
		if((Ld[i] - sum_d) < 0)
		{
			//Ld[i] = Ld[i];
			printf("LLT not done! i = %d, razn = %.15lf\n", i, Ld[i] - sum_d);
		}
		Ld[i] = sqrt(Ld[i] - sum_d);
		
	}
	
}

void MSG_my :: DI()
{
	for(int i=0; i < m; i++)
		Ld[i] = 1.0/sqrt(di[i]);

	for(int i=0; i<n_gg; i++)
		Ll[i] = 0;
}
