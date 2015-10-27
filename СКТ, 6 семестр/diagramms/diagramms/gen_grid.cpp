#define _CRT_SECURE_NO_WARNINGS 1
#include "gen_grid.h"



bool generate_2unreg_grid()
{
	int n_r, n_z;//количество подотрезков
	int m_r, m_z;//количество узлов
	double r0,rn, z0,zn, hr,hz,kr,kz;
	vector<double> r_grid_mas;
	vector<double> z_grid_mas;
	int n;
	bool flag = true;
	FILE *f = fopen("data_subseg_r.txt", "r");
	fscanf(f,"%d",&n_r);
	for(int i = 0; i < n_r; i++)
	{
		fscanf(f,"%lf %lf %lf %lf", &r0, &rn, &hr, &kr);
		if(r0 > rn ) return false; //проверка корректности данных
		
		if(kr == 1)
			n = (rn - r0)/fabs(hr) + 1;
		else n = log((rn - r0)*(kr-1)/fabs(hr) + 1)/log(kr) + 1;

		if(i == 0) flag = true;
		else flag = false;		

		generate_1unreg_grid(r0,rn,hr,kr,r_grid_mas,n,flag);
	}
	sort(r_grid_mas.begin(),r_grid_mas.end());

	f = fopen("data_subseg_z.txt", "r");
	fscanf(f,"%d",&n_z);
	for(int i = 0; i < n_z; i++)
	{
		fscanf(f,"%lf %lf %lf %lf", &z0, &zn, &hz, &kz);
		if( kz == 1)
			 n = (zn - z0)/fabs(hz) + 1;
		else  n = log((zn - z0)*(kz-1)/fabs(hz) + 1)/log(kz) + 1;

		if(i == 0) flag = true;
		else flag = false;

		generate_1unreg_grid(z0,zn,hz,kz,z_grid_mas,n, flag);
	}
	sort(z_grid_mas.begin(),z_grid_mas.end());

	m_r = r_grid_mas.size();
	m_z = z_grid_mas.size();

	printf("nodes = %d\n", m_r*m_z);
	printf("FE = %d\n", (m_r - 1)*(m_z - 1));
	f = fopen("xy.txt", "w");
	fprintf(f, "%d\n", m_r*m_z);
	for(int i=0; i<m_z; i++)
		for(int j=0; j< m_r; j++)
			fprintf(f, "%lf\t%lf\n", r_grid_mas[j], z_grid_mas[i]);
	fclose(f);

	f = fopen("nvtr_2D.txt", "w");
		fprintf(f, "%d\n", (m_r - 1)*(m_z - 1));
	for(int j = 0; j < m_z - 1; j++)
	{
		for(int i = 0; i < m_r - 1; i++)
		{
			int node_i0 = j*m_r + i;
			int node_i1 = j*m_r + i + 1;
			int node_i2 = (j+1)*m_r + i;
			int node_i3 = (j+1)*m_r + i + 1;

			fprintf(f, "%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
		}
	}
	fclose(f);

	f = fopen("ku_2D.txt", "w");
	fprintf(f, "%d\n",((m_z - 2)*2 + 2*m_r));

	for(int i=0; i<m_r; i++)
		fprintf(f,"%d\n",i);

	for(int i=1; i<m_z -1; i++)
		fprintf(f,"%d\n%d\n",i*m_r, i*m_r + (m_r -1));

	for(int i=0; i<m_r; i++)
		fprintf(f,"%d\n",m_r*(m_z -1) + i);

	fclose(f);
	return true;
}

void generate_1unreg_grid(double a, double b, double h, double kr, vector<double>&grid_mas, int n, bool flag)
{	
	double gr;
	if(a ==-3 || b==-3)
		a=a;
	if(flag)
		grid_mas.push_back(a);
	if( h > 0)
	{
		gr = a;	

		for(int i=0; i < n - 2; i++)
		{
			gr += h;
			grid_mas.push_back(gr);
			h *= kr;
		}
		gr += h;
		if(fabs((gr-b)/b) > 1E-15)
			grid_mas.push_back(gr);
		grid_mas.push_back(b);

	}
	else
	{
		gr = b;
		grid_mas.push_back(b);

		for(int i=0; i < n - 2; i++)
		{
			gr += h;
			grid_mas.push_back(gr);
			h *= kr;
		}
		gr += h;
		if(fabs((gr-a)/a) > 1E-15)
			grid_mas.push_back(gr);
	}
}


void generate_1reg_grid(double a, double b, double h, double *grid_mass, int n)
{	
	grid_mass[0] = a;
	for(int i = 1; i < n; i++)
		grid_mass[i] = grid_mass[i-1] + h;

	grid_mass[n-1] = b;
}


bool generate_2reg_grid()
{
	double a_x, a_y, b_x, b_y, h_x, h_y;
	int N, n_x, n_y;
	double *grid_mas_x, *grid_mas_y;

	FILE *f = fopen("data_reg_grid.txt", "r");
	fscanf(f,"%lf %lf %lf", &a_x, &b_x, &h_x);//, &n_x);
	fscanf(f,"%lf %lf %lf", &a_y, &b_y, &h_y);//, &n_y);
	fclose(f);

	if(a_x > b_x || h_x <= 0 || a_y > b_y || h_y <= 0) return false; //проверка корректности данных

	n_x = (b_x - a_x)/h_x + 1;
	n_y = (b_y - a_y)/h_y + 1;

	grid_mas_x = new double[n_x];
	grid_mas_y = new double[n_y];

	generate_1reg_grid(a_x, b_x, h_x, grid_mas_x, n_x);
	generate_1reg_grid(a_y, b_y, h_y, grid_mas_y, n_y);
	
	N = n_x * n_y;
	double	**grid_mas = new double* [N];
	for(int i=0; i<N; i++)
		grid_mas[i] = new double[2];

	for(int i = 0; i < n_x; i++)
	{
		for(int j = 0; j < n_y; j++)
		{
			grid_mas[j*n_x + i][0] = grid_mas_x[i];
			grid_mas[j*n_x + i][1] = grid_mas_y[j];
		}
	}
	delete[] grid_mas_x;
	delete[] grid_mas_y;

	f = fopen("xy.txt", "w");
	fprintf(f, "%d\n", N);
	for(int i=0; i<N; i++)
		fprintf(f, "%lf %lf\n", grid_mas[i][0], grid_mas[i][1]);
	fclose(f);

	f = fopen("nvtr.txt", "w");
		fprintf(f, "%d\n", (n_x - 1)*(n_y - 1));
	for(int j = 0; j < n_y - 1; j++)
	{
		for(int i = 0; i < n_x - 1; i++)
		{
			int node_i0 = j*n_x + i;
			int node_i1 = j*n_x + i + 1;
			int node_i2 = (j+1)*n_x + i;
			int node_i3 = (j+1)*n_x + i + 1;

			fprintf(f, "%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
		}
	}
	fclose(f);

	f = fopen("ku_1.txt", "w");
	fprintf(f, "%d\n",((n_y - 2)*2 + 2*n_x));
	for(int i=0; i<n_x; i++)
		fprintf(f,"%d\n",i);
	for(int i=1; i<n_y -1; i++)
		fprintf(f,"%d\n%d\n",i*n_x, i*n_x + (n_x -1));

	for(int i=0; i<n_x; i++)
		fprintf(f,"%d\n",n_x*(n_y -1) + i);

	fclose(f);

	delete[] grid_mas;
	return true;
}

