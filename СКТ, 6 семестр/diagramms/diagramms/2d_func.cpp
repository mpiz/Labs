#include "2d_func.h"


void normal_field_2D::init_iter_data(iter_data* it_data) {
	it_data->di = new double[m];
	it_data->q0 = new double[m];
	it_data->q1 = new double[m];
	it_data->q2 = new double[m];
	it_data->q_res = new double[m];
	it_data->pr = new double[m];

	it_data->gg = new double[n_gg];
	it_data->r_con = r_con;
	it_data->z_con = z0_c;

	it_data->Efi0 = new double[m];
	it_data->Efi1 = new double[m];
	it_data->Efi2 = new double[m];
	it_data->Efi_res = new double[m];

}

void normal_field_2D::read_ku()
{
	FILE *f;
	f = fopen("ku_2D.txt", "r");
	fscanf(f, "%d", &k1);
	int tmp;
	ku_one = new ku_o[k1];
	for (int i = 0; i<k1; i++)
	{
		fscanf(f, "%d", &tmp);
		onekuset.insert(tmp);
		ku_one[i].node = tmp;
	}
	fclose(f);
}

void normal_field_2D::read_KE()
{
	FILE *f = fopen("xy.txt", "r");
	fscanf(f, "%d ", &m);
	nodes = new node[m];
	for (int i = 0; i<m; i++)
		fscanf(f, "%lf %lf", &nodes[i].r, &nodes[i].z);

	fclose(f);

	f = fopen("nvtr_2D.txt", "r");
	fscanf(f, "%d ", &N);
	el = new KE[N];
	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<4; j++)
			fscanf(f, "%d ", &el[i].node[j]);
	}
	fclose(f);

	//form_file_nkvat();
}
void normal_field_2D::form_file_nkvat()
{
	FILE *f = fopen("nvkat.txt", "w");
	FILE *f1 = fopen("nvkats.txt", "w");
	fprintf(f, "%d\n", N);
	fprintf(f1, "%d\n", N);
	for (int i = 0; i<N; i++)
	{
		el[i].area = form_area(nodes[el[i].node[0]].r, nodes[el[i].node[0]].z, nodes[el[i].node[2]].z);
		fprintf(f, "%d\n", el[i].area);
		fprintf(f1, "%d\n", 1);
	}
	fclose(f);
	fclose(f1);


}

void normal_field_2D::read_time()
{
	FILE *f = fopen("data_time.txt", "r");
	double h_t, t0, tn, k_r;
	fscanf(f, "%lf %lf", &t0, &tn);
	fscanf(f, "%lf", &h_t);
	fscanf(f, "%lf", &k_r);

	if (k_r == 1)
		n_time = (int)((tn - t0) / fabs(h_t) + 1);
	else n_time = log((tn - t0)*(k_r - 1) / fabs(h_t) + 1) / log(k_r) + 1;


	time_mas = new double[n_time];

	time_mas[0] = t0;
	for (int i = 1; i < n_time; i++)
	{
		time_mas[i] = time_mas[i - 1] + h_t;
		h_t *= k_r;
	}
	time_mas[n_time - 1] = tn;
	fclose(f);
	f = fopen("time_mas.txt", "w");
	fprintf(f, "%d\n", n_time);
	for (int i = 0; i<n_time; i++)
		fprintf(f, "%.8e\n", time_mas[i]);
	fclose(f);
}

void normal_field_2D::read()
{
	FILE *f;
	read_KE();
	read_ku();
	read_time();
	//выделение памяти под матрицу m - размерность матрицы
	f = fopen("portret_2D.txt", "r");
	fscanf(f, "%d", &n_gg);
	ig = new int[m + 1];
	jg = new int[n_gg];

	for (int i = 0; i<m + 1; i++)
		fscanf(f, "%d", &ig[i]);

	for (int i = 0; i<n_gg; i++)
		fscanf(f, "%d", &jg[i]);
	fclose(f);

	read_move_controle();
}

double normal_field_2D::form_F(double r, double z, int i, int t_i)
{
	double f;
	double t = time_mas[t_i];
	double sigma = form_sigma(r, z, i, t_i);
	/*if(nodes[i].r == 0.5 && nodes[i].z == 0)
	{
	f = 10;
	printf("%d\n", i);
	}
	else f = 0;*/
	//f = t*t/(r*r) + sigma*2*t;
	//f = 2*r*z*t;
	//f = r*z;
	f = 0;
	//double mu = form_mu(r,z,i,t_i);
	//f = -3*mu;
	//	f = -3*t*mu;
	//f = 3.0*r*z*t*t;
	//f = r*z;
	return f;
}

double normal_field_2D::form_u(double r, double z, int i, int t_i)
{
	double t = time_mas[t_i];
	//double u = t*t;
	//double u = r*z*t;
	//double u = r*r;
	double u = r*z*t;
	//double u = r*z;
	return u;
}

double normal_field_2D::form_sigma(double r, double z, int i, int t_i)
{		
		if ((z > 0)) return 0;
		if ((z <= 0 && z > sigma_sol[0])) return sigma_sol[1];
		if ((z <= sigma_sol[0] && z > sigma_sol[2])) return sigma_sol[3];
		if (z <= sigma_sol[2] && z > sigma_sol[4]) return sigma_sol[5];
		if (z <= sigma_sol[4] && z > sigma_sol[6]) return sigma_sol[7];
		if (z <= sigma_sol[7]) return 0.1;
	
}


int normal_field_2D::form_area(double r, double z1, double z2)
{
//	int area;
	if (r >= 1.1)
	{
		if ((z1 + z2) / 2 > 0) return 1;
		if ((z1 + z2) / 2 < -6 && (z1 + z2) / 2 > -11) return 2;
		if ((z1 + z2) / 2 <= 11 && (z1 + z2) / 2 > -19) return 3;
		if ((z1 + z2) / 2 <= -19 && (z1 + z2) / 2 > -23) return 6;
		return 4;
	}
	else return 5;
	

}

double normal_field_2D::form_mu(double r, double z, int i, int t_i)
{
	double mu = 100;
	double mu0 = 16 * atan(1.0)*1e-7;
	return 1.0 / (mu0);

}
void normal_field_2D::form_local_matrix_3sl(int i, int t_i, double L[4][4], double b[4], iter_data* it_data)
{

	double g = 0;
	double p[4];
	double dt0 = time_mas[t_i] - time_mas[t_i - 1];
	double dt1 = time_mas[t_i - 1] - time_mas[t_i - 2];
	double dt = time_mas[t_i] - time_mas[t_i - 2];
	for (int j = 0; j<4; j++)
	{
		p[j] = form_F(nodes[el[i].node[j]].r, nodes[el[i].node[j]].z, i, t_i);
		g += form_sigma(nodes[el[i].node[j]].r, nodes[el[i].node[j]].z, i, t_i);
	}
	g = g / 4.0;
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;

	//матрица массы и жесткости
	double M[4][4];
	double G[4][4];
	double D[4][4];
	//form_matrix_M(nodes[el[i].node[0]].r,M,hr,hz);
	form_matrix_M1(nodes[el[i].node[0]].r, M, hr, hz);

	//матрица жесткости
	//form_matrix_G(hr,hz,G,i,nodes[el[i].node[0]].r,t_i);
	//form_matrix_G1(hr, hz, G, nodes[el[i].node[0]].r);
	form_matrix_G2(hr, hz, nodes[el[i].node[0]].r, nodes[el[i].node[1]].r, G);

	//матрица добавки
	//form_matrix_D(nodes[el[i].node[0]].r,hr,hz,D);
	form_matrix_D1(nodes[el[i].node[0]].r, hr, hz, D);

	double mu = form_mu(nodes[el[i].node[0]].r, nodes[el[i].node[0]].z, i, t_i);
	//формируем локальную матрицу и правую часть
	for (int j = 0; j<4; j++)
	{
		b[j] = 0;
		for (int k = 0; k<4; k++)
		{
			L[j][k] = g*M[j][k] * (dt + dt0) / (dt*dt0) + mu*G[j][k] + mu*D[j][k];
			b[j] += g*M[j][k] * it_data->q0[el[i].node[k]] * (-dt0 / (dt*dt1)) + g*M[j][k] * it_data->q1[el[i].node[k]] * (dt / (dt0*dt1)) + M[j][k] * p[k];

		}


	}
}
void normal_field_2D::form_local_matrix_2sl(int i, int t_i, double L[4][4], double b[4], iter_data* it_data)
{
	double g = 0;
	double p[4];
	double dt = time_mas[t_i] - time_mas[t_i - 1];
	for (int j = 0; j<4; j++)
	{
		p[j] = form_F(nodes[el[i].node[j]].r, nodes[el[i].node[j]].z, i, t_i);
		g += form_sigma(nodes[el[i].node[j]].r, nodes[el[i].node[j]].z, i, t_i);
	}
	g = g / 4.0;
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;

	//матрица массы и жесткости
	double M[4][4];
	double G[4][4];
	double D[4][4];

	form_matrix_M1(nodes[el[i].node[0]].r, M, hr, hz);

	//матрица жесткости
	form_matrix_G2(hr, hz, nodes[el[i].node[0]].r, nodes[el[i].node[1]].r, G);

	//матрица добавки
	form_matrix_D1(nodes[el[i].node[0]].r, hr, hz, D);
	double mu = form_mu(nodes[el[i].node[0]].r, nodes[el[i].node[0]].z, i, t_i);

	//формируем локальную матрицу и правую часть
	for (int j = 0; j<4; j++)
	{
		b[j] = 0;
		for (int k = 0; k<4; k++)
		{
			L[j][k] = g*M[j][k] / dt + mu*G[j][k] + mu*D[j][k];
			b[j] += g*M[j][k] * it_data->q0[el[i].node[k]] / dt + M[j][k] * p[k];
		}
	}

}
void normal_field_2D::form_local_matrix_stat(int i, int t_i, double L[4][4], double b[4])
{

	double p[4];
	for (int j = 0; j<4; j++)
		p[j] = form_F(nodes[el[i].node[j]].r, nodes[el[i].node[j]].z, i, t_i);


	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;

	//матрица массы и жесткости
	double M[4][4];
	double G[4][4];
	double D[4][4];

	//матрица жесткости
	//form_matrix_G(hr,hz,G,i,nodes[el[i].node[0]].r,t_i);
	//form_matrix_G1(hr, hz, G, nodes[el[i].node[0]].r);
	form_matrix_G2(hr, hz, nodes[el[i].node[0]].r, nodes[el[i].node[1]].r, G);

	//матрица массы
	//form_matrix_M(nodes[el[i].node[0]].r,M,hr,hz);
	form_matrix_M1(nodes[el[i].node[0]].r, M, hr, hz);

	//матрица добавки
	//form_matrix_D(nodes[el[i].node[0]].r,hr,hz,D);
	form_matrix_D1(nodes[el[i].node[0]].r, hr, hz, D);

	double mu = form_mu(nodes[el[i].node[0]].r, nodes[el[i].node[0]].z, i, t_i);

	//формируем локальную матрицу и правую часть
	for (int j = 0; j<4; j++)
	{
		b[j] = 0;
		for (int k = 0; k<4; k++)
		{
			L[j][k] = mu*G[j][k] + mu*D[j][k];
			b[j] += M[j][k] * p[k];

		}
	}

}

void normal_field_2D::form_global_matrix(int t_i, iter_data* it_data)
{
	int k;
	double L[4][4], b[4];

	//матрица без краевых условий
	for (int kon = 0; kon<N; kon++)
	{
		if (t_i == 0) form_local_matrix_stat(kon, t_i, L, b);
		if (t_i == 1) form_local_matrix_2sl(kon, t_i, L, b, it_data);
		if (t_i > 1) form_local_matrix_3sl(kon, t_i, L, b, it_data);
		for (int i = 0; i<4; i++)
		{
			for (int j = 0; j<i; j++)
			{
				int nach = ig[el[kon].node[i]];
				//	int kol = ig[el[kon].node[i] + 1] - nach;
				bool not_end = true;

				for (k = nach; not_end; k++)
					if (jg[k] == el[kon].node[j]) not_end = false;

				it_data->gg[k - 1] += L[i][j];
			}
			it_data->di[el[kon].node[i]] += L[i][i];
			it_data->pr[el[kon].node[i]] += b[i];

		}
	}


	if (t_i == 0) 
		controle(it_data);
	ku(t_i, it_data);
	

}
/*
void normal_field_2D::mult_matrix(double *x, double *&y)
{
	for (int i = 0; i<m; i++)
	{
		y[i] = di[i] * x[i];

		for (int j = ig[i]; j<ig[i + 1]; j++)
		{
			y[jg[j]] += gg[j] * x[i];
			y[i] += gg[j] * x[jg[j]];
		}
	}
}
*/

void normal_field_2D::ku(int t_i, iter_data* it_data)
{
	int k;
	for (k = 0; k<m; k++)
	{// пробегаем по всей матрице
		if (onekuset.find(k) != onekuset.end())
		{
			//double val = form_u(nodes[k].r, nodes[k].z, 0, t_i);;//значение в узле
			double val = 0;
			it_data->pr[k] = val;
			it_data->di[k] = 1;

			for (int i = ig[k]; i<ig[k + 1]; i++)
			{
				if (onekuset.find(jg[i]) == onekuset.end())
					it_data->pr[jg[i]] -= it_data->gg[i] * val;
					it_data->gg[i] = 0;
			}
		}
		else
		{
			for (int i = ig[k]; i<ig[k + 1]; i++)
			{
				if (onekuset.find(jg[i]) != onekuset.end())
				{
					//double val = form_u(nodes[jg[i]].r, nodes[jg[i]].z, 0 ,t_i);//значение в узле
					double val = 0;
					it_data->pr[k] -= it_data->gg[i] * val;
					it_data->gg[i] = 0;
				}
			}
		}
	}
	
}

void normal_field_2D::controle(iter_data* it_data)
{
	double u = r_con;
	bool not_end = true;

	for (int i = 0; i<N && not_end; i++)
	{
		if ((nodes[el[i].node[0]].r <= r_con) &&
			(nodes[el[i].node[1]].r >= r_con) &&
			(nodes[el[i].node[0]].z <= it_data->z_con) &&
			(nodes[el[i].node[2]].z >= it_data->z_con))
		{


			it_data->pr[el[i].node[0]] += u*fi_1(i, r_con, it_data->z_con);
			it_data->pr[el[i].node[1]] += u*fi_2(i, r_con, it_data->z_con);
			it_data->pr[el[i].node[2]] += u*fi_3(i, r_con, it_data->z_con);
			it_data->pr[el[i].node[3]] += u*fi_4(i, r_con, it_data->z_con);
			not_end = false;
		}
	}
}

void normal_field_2D::calc_u(int t_i, FILE *fi, FILE *f1, double *qq)
{

	double norma = 0;
	double norma_u = 0;
	u = new double[m];
	//fprintf(fi, "time =%d\n", t_i);
	fflush(fi);
	//fprintf(f1, "time =%d\n", t_i);

	for (int i = 0; i<m; i++)
	{
		//	int i = ku_one[j].node;
		u[i] = form_u(nodes[i].r, nodes[i].z, 0, t_i);
		fprintf(fi, "%d\t%.15lf\t%.15lf\n", i, u[i], qq[i]);
		fflush(fi);
		norma_u += u[i] * u[i];
		u[i] -= qq[i];
		norma += u[i] * u[i];

	}
	fflush(fi);
	norma = sqrt(norma / norma_u);
	//fprintf(fi, "%d\ttime=\t%lf\tnorma=\t%.8e\n",t_i, time_mas[t_i], norma);
	//if(t_i % 4 == 0)
	fprintf(fi, "%.8e\t%.8e\n", time_mas[t_i], norma);
	fflush(fi);
}

void normal_field_2D::printf_matrix(int i)
{
	FILE *f, *ff, *fff, *ffff;
	char name[40];
	sprintf(name, "matrix%d.txt", i);
	f = fopen(name, "w");

	sprintf(name, "prch%d.txt", i);
	ff = fopen(name, "w");

	sprintf(name, "di%d.txt", i);
	fff = fopen(name, "w");

	sprintf(name, "q%d.txt", i);
	ffff = fopen(name, "w");

	fprintf(f, "%d\n", n_gg);
	fprintf(ff, "%d\n", m);
	fprintf(fff, "%d\n", m);
	fprintf(ffff, "%d\n", m);


	fclose(f);
	fclose(ff);
	fclose(fff);
	fclose(ffff);

}

void normal_field_2D::move_inverse()
{
	double *b1, **A, ***der, **dif, *dsigma, **A_tmp;
	double *find_sigma = new double[2*n_sigma];

	double *find_sigma_cond = new double[2 * n_sigma];

	FILE *f = fopen("first_sigma.txt", "r");

	for (int i = 0; i < 2*n_sigma; i++)
		fscanf(f, "%lf", &find_sigma[i]);
	fclose(f);
	double F_new = 0, F_old = 0;
	double **m_rec1 = new double *[n_time];
	for (int i = 0; i < n_time; i++) {
		m_rec1[i] = new double[n_con];
		for (int j = 0; j < n_con; j++)
			m_rec1[i][j] = 0;
	}

	moving_controle(find_sigma, m_rec1);//считаем прямую задачу в них же сразу находим значения rec1
	//и находим значение функционала

	for (int i = 0; i < n_time; i++)
		for (int j = 0; j < n_con; j++)
			F_new += (m_rec_true[i][j] - m_rec1[i][j])*(m_rec_true[i][j] - m_rec1[i][j]) / (m_rec_true[i][j] * m_rec_true[i][j]);


	b1 = new double[2 * n_sigma];
	A = new double*[2 * n_sigma];
	A_tmp = new double*[2 * n_sigma];
	der = new double**[2 * n_sigma];//i - параметр

	for (int i = 0; i < 2 * n_sigma; i++)
	{
		A[i] = new double[2 * n_sigma];
		A_tmp[i] = new double[2 * n_sigma];
		der[i] = new double*[n_time];//j - время

		for (int j = 0; j < n_time; j++)
			der[i][j] = new double[n_con];//k - источник-приемник
	}

	dif = new double*[n_time];
	for (int i = 0; i < n_time; i++)
		dif[i] = new double[n_con];

	dsigma = new double[2 * n_sigma];

	bool not_end = true;

	if (F_new > 1E-15)
	{
		not_end = false;
	}


	//rec2 для приращения, rec1 текущая сигма
	while (!not_end)
	{
		F_old = F_new;

		//cобираем матрицу
#pragma omp parallel for
		for (int i = 0; i < 2 * n_sigma; i++)
		{
			double** m_rec2 = new double*[n_time];
			double *find_sigma1 = new double[2 * n_sigma];

			for (int i = 0; i < n_time; i++)
					m_rec2[i] = new double[n_con];	
			
			for (int pp = 0; pp < 2 * n_sigma; pp++)
				find_sigma1[pp] = find_sigma[pp];

			double ds = 0.1*find_sigma[i];
			find_sigma1[i] += ds;

			for (int i = 0; i < n_time; i++)
				for (int j = 0; j < n_con; j++)
					m_rec2[i][j] = 0;

			while ((moving_controle(find_sigma1, m_rec2) == 0))
			{
				find_sigma1[i] -= 0.5*ds;
			}

			printf("rec2");

			for (int j = 0; j < n_time; j++)
				for (int k = 0; k < n_con; k++)
				{
				double val = (m_rec1[j][k] - m_rec2[j][k]) / ds;
				der[i][j][k] = val;
				}


			for (int i = 0; i < n_time; i++)
				delete m_rec2[i];

			delete m_rec2;
			delete find_sigma1;
		}

		for (int j = 0; j < n_time; j++)
			for (int k = 0; k < n_con; k++)
				dif[j][k] = m_rec_true[j][k] - m_rec1[j][k];


		for (int i = 0; i < 2 * n_sigma; i++)
		{
			b1[i] = 0;
			for (int j = 0; j < n_time; j++)
				for (int k = 0; k < n_con; k++)
					b1[i] -= (der[i][j][k] * dif[j][k]) / (m_rec_true[j][k] * m_rec_true[j][k]);

			for (int i1 = 0; i1 < 2 * n_sigma; i1++)
			{
				A[i][i1] = 0;
				A_tmp[i][i1] = 0;
				for (int j = 0; j < n_time; j++)
					for (int k = 0; k < n_con; k++)
						A[i][i1] += (der[i][j][k] * der[i1][j][k]) / (m_rec_true[j][k] * m_rec_true[j][k]);

				A_tmp[i][i1] = A[i][i1];
			}
		}


		// тут надо решить слау A*dsigma = b
		double *alfa = new double[2 * n_sigma];
		for (int i = 0; i < 2 * n_sigma; i++)
		{
			alfa[i] = 1E-15;
			A_tmp[i][i] += alfa[i];
		}
		SLAE_solution_Gauss(A_tmp, b1, dsigma, 2 * n_sigma);

		for (int i1 = 0; i1 < 2 * n_sigma; i1++)
			for (int j1 = 0; j1 < 2 * n_sigma; j1++)
				A_tmp[i1][j1] = A[i1][j1];

		bool alfa_end = false;

		while (!alfa_end)
		{
			int kol = 0;
			for (int i = 1; i < 2 * n_sigma; i = i + 2)
			{
				if (!(0 < (find_sigma[i] + dsigma[i]) && find_sigma[i] + dsigma[i] < 10))
				{
					alfa[i] = alfa[i] * 2;
					A_tmp[i][i] += alfa[i];
					kol++;
				}

				if (alfa[i] > 1e-10)
				{
					alfa_end = true;

				}
			} 

			for (int i = 2; i < 2 * n_sigma; i = i + 2)
			{
				if ((find_sigma[i - 2] < (find_sigma[i] + dsigma[i])))
				{
					alfa[i] = alfa[i] * 2;
					A_tmp[i][i] += alfa[i];
					kol++;
				}

				if (alfa[i] > 1e-10)
				{
					alfa_end = true;
				}
				
			}

			if (kol == 0)
			{
				alfa_end = true;
				//SLAE_solution_Gauss(A, b1, dsigma, kol_param);
			}
			else
			{
				if (alfa_end != true)
				{
					SLAE_solution_Gauss(A_tmp, b1, dsigma, 2 * n_sigma);
					for (int i1 = 0; i1 < 2 * n_sigma; i1++)
						for (int j1 = 0; j1 < 2 * n_sigma; j1++)
							A_tmp[i1][j1] = A[i1][j1];
				}


			}

		}


		for (int i = 0; i < 2 * n_sigma; i++)
		{
			if (i%2 == 0) dsigma[i] = int(dsigma[i] * 10.0) / 10.0;
			find_sigma_cond[i] = find_sigma[i] + dsigma[i];
		}


		bool cont = true;
		do {

			for (int i = 0; i < n_time; i++)
				for (int j = 0; j < n_con; j++)
					m_rec1[i][j] = 0;
			if (moving_controle(find_sigma_cond, m_rec1) == 1)//если все хорошо
			{
				F_new = 0;
				for (int i = 0; i < n_time; i++)
					for (int j = 0; j < n_con; j++)
						F_new += (m_rec_true[i][j] - m_rec1[i][j])*(m_rec_true[i][j] - m_rec1[i][j]) / (m_rec_true[i][j] * m_rec_true[i][j]);

				if (F_new <= F_old)
					cont = false;
				else
				{
					for (int i = 0; i < 2 * n_sigma; i++)
					{
						dsigma[i] = 0.5*dsigma[i];
						find_sigma_cond[i] = find_sigma[i] + dsigma[i];
					}
				}
			}
			else
			{

				for (int i = 0; i < 2 * n_sigma; i++)
				{
					dsigma[i] = 0.5*dsigma[i];
					find_sigma_cond[i] = find_sigma[i] + dsigma[i];
				}
			}

		} while (cont);


		FILE *f_log = fopen("f_log.txt", "a");
		for (int i = 0; i < 2 * n_sigma; i++)
			find_sigma[i] = find_sigma_cond[i];
		fprintf(f_log, "F:\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n", F_old, F_new, dsigma[0], dsigma[1], dsigma[2]);

		for (int i = 0; i < 2 * n_sigma; i++)
			fprintf(f_log, "sigma%d = %.15lf ", i, find_sigma[i]);

		fprintf(f_log, "\n");
		fclose(f_log);

		//printf("F:\t%.3e\t%.3e\t%.3e\t%.3et%.3e\n", F_old, F_new, dsigma[0], dsigma[1], dsigma[2]);// , dsigma[3]);
		//printf("%.15lf\t%.15lf\t%.15lf\n", find_sigma[number_param[0]], find_sigma[number_param[1]], find_sigma[number_param[2]]);//, find_sigma[number_param[2]], find_sigma[number_param[3]]);
		if ((F_new < 1E-15) || (fabs(F_new - F_old) < 1E-15))
		{
			not_end = true;

		}

	}

	for (int i = 0; i < n_con; i++)
		delete m_rec1[i];

	delete m_rec1;
	delete find_sigma_cond;
}

void normal_field_2D::move_data_true()
{

	FILE *f;

	//считываем истинную сигму и считаем значения в приемниках это данные
	f = fopen("true_sigma.txt", "r");
	fscanf(f, "%d", &n_sigma);
	true_sigma = new double[2*n_sigma];
	sigma_sol = new double[2*n_sigma];

	for (int i = 0; i<2*n_sigma; i++)
		fscanf(f, "%lf", &true_sigma[i]);
	fclose(f);

	FILE *f_time = fopen("time_mas.txt", "r");
	fscanf(f_time, "%d", &n_time);
	fclose(f_time);

	m_rec_true = new double*[n_time];
	
	for (int j = 0; j < n_time; j++) {
		m_rec_true[j] = new double[n_con];
		for (int k = 0; k < n_con; k++)
			m_rec_true[j][k] = 0;
	}
	moving_controle(true_sigma, m_rec_true);

#ifdef NOISE

	double r1 = rand() * 1.0 / RAND_MAX;
	double r2 = rand() * 1.0 / RAND_MAX;

	double pi = 4 * atan(1.0);

	double norm = cos(2 * pi*r1) * sqrt(-2 * log(r2));
	double w = 0.05;

	
	for (int j = 0; j < n_time; j++)
		for (int k = 0; k < n_rec; k++)
		{
			if (k == 0)
			m_rec_true[j][k] = m_rec_true[j][k] + w*norm * m_rec_true[j][k];
		}

#endif
}


void normal_field_2D::time_calc_matrix(int i_con, double *fsigma, double ** rec, iter_data* it_data)
{
//	FILE *ff;
	printf("nu begin\n");
	int th_i = omp_get_thread_num();
	
	MSG_my s_MSG;


	for (int i = 0; i<n_time; i++)
	{

		for (int j = 0; j<n_gg; j++)
			it_data->gg[j] = 0;
		for (int j = 0; j<m; j++)
		{
			it_data->di[j] = 0;
			it_data->pr[j] = 0;
		}

		if (i == 0)
		{

			
			for (int j = 0; j<m; j++)
				it_data->q0[j] = 0.0;

			form_global_matrix(i, it_data);			
			s_MSG.give_data(ig, jg, it_data->gg, it_data->di, it_data->pr, it_data->q0, m, n_gg);
			Efi(i, it_data);
			for (int j = 0; j<m; j++)
			{
				it_data->Efi_res[j] = it_data->Efi0[j];
				it_data->q_res[j] = it_data->q0[j];
			}
			double val = Efi_inpoint(r_con, it_data->z_con - 4, it_data);
			rec[i][i_con] = val;
			printf("%d 2D calc q0\n", th_i);
			


		}
		if (i > 1)
		{
					
			for (int j = 0; j<m; j++)
			{
				it_data->q2[j] = 0.0;
				//x0[j] = 0.0;
			}
			form_global_matrix(i, it_data);

			s_MSG.give_data(ig, jg, it_data->gg, it_data->di, it_data->pr, it_data->q2, m, n_gg);
			Efi(i, it_data);
			for (int j = 0; j<m; j++)
			{
				it_data->Efi_res[j] = it_data->Efi2[j];
				it_data->q_res[j] = it_data->q2[j];
			}
			double val = Efi_inpoint(r_con, it_data->z_con - 4, it_data);
			rec[i][i_con] = val;
			for (int ii = 0; ii<m; ii++)
			{
				it_data->q0[ii] = it_data->q1[ii];
				it_data->q1[ii] = it_data->q2[ii];
				it_data->q2[ii] = 0;
			}

			if (i % 30 == 0 || i == n_time-1)
				printf("%d 2D calc q%d\n", th_i, i);
			


		}
		if (i == 1)
		{
			for (int j = 0; j<m; j++)
			{
				it_data->q1[j] = 0;
				//x0[j] = 0;
			}
			form_global_matrix(i, it_data);

			s_MSG.give_data(ig, jg, it_data->gg, it_data->di, it_data->pr, it_data->q1, m, n_gg);
			Efi(i, it_data);

			for (int j = 0; j<m; j++)
			{
				it_data->Efi_res[j] = it_data->Efi1[j];
				it_data->q_res[j] = it_data->q1[j];
			}
			double val = Efi_inpoint(r_con, it_data->z_con - 4, it_data);
			rec[i][i_con] = val;
			printf("%d 2D calc q1\n", th_i);

		}
	}
	
	printf("nu end\n");
}

int normal_field_2D::moving_controle(double *fsigma, double ** rec)
{
	for (int i = 0; i < 2 * n_sigma; i++)
		sigma_sol[i] = fsigma[i];

	iter_data it_data;
	init_iter_data(&it_data);

	int th_i = omp_get_thread_num();

	for (int i = 0; i < n_con; i++)
	{ 
		printf("%d started calc z = %lf\n", th_i, it_data.z_con);
		time_calc_matrix(i, fsigma, rec, &it_data);
		it_data.z_con -= hcon;
	}

	return 1;
}

void normal_field_2D :: read_move_controle()
{
	FILE *f = fopen("data_controle.txt", "r");
	
	fscanf(f, "%lf %lf %lf %d", &z0_c, &r_con, &hcon, &n_con);
	fclose(f);

}

double normal_field_2D::fi_1(int i, double r, double z)
{
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;
	double fi_1 = ((nodes[el[i].node[1]].r - r) / hr) * ((nodes[el[i].node[2]].z - z) / hz);
	return fi_1;
}
double normal_field_2D::fi_2(int i, double r, double z)
{
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;
	double fi_2 = ((r - nodes[el[i].node[0]].r) / hr) * ((nodes[el[i].node[2]].z - z) / hz);
	return fi_2;
}
double normal_field_2D::fi_3(int i, double r, double z)
{
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;
	double fi_3 = ((nodes[el[i].node[1]].r - r) / hr) * ((z - nodes[el[i].node[0]].z) / hz);
	return fi_3;
}
double normal_field_2D::fi_4(int i, double r, double z)
{
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;
	double fi_4 = ((r - nodes[el[i].node[0]].r) / hr) * ((z - nodes[el[i].node[0]].z) / hz);
	return fi_4;
}

double normal_field_2D::pr_fi_1(int i, double r, double z)
{
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;
	double fi_1 = (-1.0) / hr * (nodes[el[i].node[2]].z - z) / hz;
	return fi_1;
}
double normal_field_2D::pr_fi_2(int i, double r, double z)
{
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;
	double fi_2 = 1.0 / hr * (nodes[el[i].node[2]].z - z) / hz;
	return fi_2;
}
double normal_field_2D::pr_fi_3(int i, double r, double z)
{
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;
	double fi_3 = (-1.0) / hr * (z - nodes[el[i].node[0]].z) / hz;
	return fi_3;
}
double normal_field_2D::pr_fi_4(int i, double r, double z)
{
	double hr = nodes[el[i].node[1]].r - nodes[el[i].node[0]].r;
	double hz = nodes[el[i].node[2]].z - nodes[el[i].node[0]].z;
	double fi_4 = 1.0 / hr * (z - nodes[el[i].node[0]].z) / hz;
	return fi_4;
}

double normal_field_2D::Afi_inpoint(double r, double z, iter_data* it_data)
{
	int i = 0;
	bool not_end = true;

	for (i = 0; i<N && not_end; i++)
	{
		if ((nodes[el[i].node[0]].r <= r) &&
			(nodes[el[i].node[1]].r >= r) &&
			(nodes[el[i].node[0]].z <= z) &&
			(nodes[el[i].node[2]].z >= z)) not_end = false;
	}
	i--;

	double Afi = (it_data->q_res[el[i].node[0]] * fi_1(i, r, z) + it_data->q_res[el[i].node[1]] * fi_2(i, r, z) + it_data->q_res[el[i].node[2]] * fi_3(i, r, z) + it_data->q_res[el[i].node[3]] * fi_4(i, r, z));

	return Afi;
}

int normal_field_2D::find_KE(double r, double z)
{
	int i = 0;
	bool not_end = true;

	for (i = 0; i<N && not_end; i++)
	{
		if ((nodes[el[i].node[0]].r <= r) &&
			(nodes[el[i].node[1]].r >= r) &&
			(nodes[el[i].node[0]].z <= z) &&
			(nodes[el[i].node[2]].z >= z)) not_end = false;
	}
	i--;

	return i;
}

double normal_field_2D::Efi_inpoint(double r, double z, iter_data* it_data)
{
	int i = find_KE(r, z);

	double Efi = (it_data->Efi_res[el[i].node[0]] * fi_1(i, r, z) + it_data->Efi_res[el[i].node[1]] * fi_2(i, r, z) + it_data->Efi_res[el[i].node[2]] * fi_3(i, r, z) + it_data->Efi_res[el[i].node[3]] * fi_4(i, r, z));

	return Efi;
}


void normal_field_2D::Efi(int t_i, iter_data* it_data)
{
	
		if (t_i == 0)
		{
			for (int j = 0; j<m; j++)
				it_data->Efi0[j] = it_data->q0[j];
		}
		if (t_i == 1)
		{
			double dt = time_mas[t_i] - time_mas[t_i - 1];
			for (int j = 0; j < m; j++)
				it_data->Efi1[j] = -(it_data->q1[j] / dt - it_data->q0[j] / dt);

		}
		if (t_i > 1)
		{
			for (int j = 0; j<m; j++)
				it_data->Efi2[j] = 0;
			double dt0 = time_mas[t_i] - time_mas[t_i - 1];
			double dt1 = time_mas[t_i - 1] - time_mas[t_i - 2];
			double dt = time_mas[t_i] - time_mas[t_i - 2];

			for (int j = 0; j<m; j++)
				it_data->Efi2[j] = -((dt0 / (dt1*dt))* it_data->q0[j] - (dt / (dt0*dt1))* it_data->q1[j] + ((dt + dt0) / (dt*dt0))* it_data->q2[j]);
			
			double *c;
			
			
		}

}



//сборка матриц М и D G


int normal_field_2D::mulin(int i)
{
	return i % 2;
}
int normal_field_2D::nulin(int j)
{
	return j / 2;
}
void normal_field_2D::form_matrix_D1(double r, double hr, double hz, double D[4][4])
{
	double d = r / hr;
	double D0[2][2];
	double M0[2][2];
	double ln = log(1 + 1.0 / d);
	D0[0][0] = ln * (1 + d) * (1 + d) - d - 1.5;
	D0[0][1] = -ln * (1 + d) * d + d + 0.5;
	D0[1][0] = -ln * (1 + d) * d + d + 0.5;
	D0[1][1] = ln * d * d - d + 0.5;


	double C = hz / 6.0;
	M0[0][0] = 2.0 * C;
	M0[0][1] = C;
	M0[1][0] = C;
	M0[1][1] = 2.0 * C;

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			D[i][j] = M0[nulin(i)][nulin(j)] * D0[mulin(i)][mulin(j)];
}

void normal_field_2D::form_matrix_M1(double rp, double M[4][4], double hr, double hz)
{
	double C = hz / 6.0;
	double M0[2][2];
	M0[0][0] = 2.0 * C;
	M0[0][1] = C;
	M0[1][0] = C;
	M0[1][1] = 2.0 * C;

	double c1 = hr / 6.0;
	double c2 = hr / 2.0;
	double R0[2][2];
	R0[0][0] = c1 * (2 * rp + c2);
	R0[0][1] = c1 * (rp + c2);
	R0[1][0] = c1 * (rp + c2);
	R0[1][1] = c1 * (2 * rp + 3 * c2);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			M[i][j] = M0[nulin(i)][nulin(j)] * R0[mulin(i)][mulin(j)];
}


void normal_field_2D::form_matrix_G2(double hr, double hz, double r1, double r2, double G[4][4])
{
	double c1 = hr / 6.0;
	double c2 = hr / 2.0;
	double R0[2][2];
	R0[0][0] = c1 * (2 * r1 + c2);
	R0[0][1] = c1 * (r1 + c2);
	R0[1][0] = c1 * (r1 + c2);
	R0[1][1] = c1 * (2 * r1 + 3 * c2);

	double C = hz / 6.0;
	double M0[2][2];
	M0[0][0] = 2.0 * C;
	M0[0][1] = C;
	M0[1][0] = C;
	M0[1][1] = 2.0 * C;

	C = 1.0 / hz;
	double GZ[2][2];
	GZ[0][0] = C;
	GZ[0][1] = -C;
	GZ[1][0] = -C;
	GZ[1][1] = C;

	double GR[2][2];
	GR[0][0] = (r2 + r1) / (2 * (r2 - r1));
	GR[0][1] = -(r2 + r1) / (2 * (r2 - r1));
	GR[1][0] = -(r2 + r1) / (2 * (r2 - r1));
	GR[1][1] = (r2 + r1) / (2 * (r2 - r1));

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
			G[i][j] = (M0[nulin(i)][nulin(j)] * GR[mulin(i)][mulin(j)] + GZ[nulin(i)][nulin(j)] * R0[mulin(i)][mulin(j)]);
}
