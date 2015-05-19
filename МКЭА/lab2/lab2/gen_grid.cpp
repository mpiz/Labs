#include "gen_grid.h"

bool generate_3unreg_grid()
{
	int n_x, n_y, n_z;//количество подотрезков
	int m_x, m_y, m_z;//количество узлов
	double x0,xn,y0,yn, z0,zn, hx, hy, hz,kx, ky,kz;
	vector<double> x_grid_mas;
	vector<double> z_grid_mas;
	vector<double> y_grid_mas;
	int n;
	bool flag = true;
	FILE *f = fopen("data_subseg_x.txt", "r");
	fscanf(f,"%d",&n_x);
	for(int i = 0; i < n_x; i++)
	{
		fscanf(f,"%lf %lf %lf %lf", &x0, &xn, &hx, &kx);
		if(x0 > xn ) return false; //проверка корректности данных
		
		if(kx == 1)
			n = (xn - x0)/fabs(hx) + 1;
		else n = log((xn - x0)*(kx-1)/fabs(hx) + 1)/log(kx) + 1;

		if(i == 0) flag = true;
		else flag = false;		

		generate_1unreg_grid(x0,xn,hx,kx,x_grid_mas,n,flag);
	}
	sort(x_grid_mas.begin(),x_grid_mas.end());

	f = fopen("data_subseg_z3D.txt", "r");
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
	

	f = fopen("data_subseg_y.txt", "r");
	fscanf(f,"%d",&n_y);
	for(int i = 0; i < n_y; i++)
	{
		fscanf(f,"%lf %lf %lf %lf", &y0, &yn, &hy, &ky);
		if( ky == 1)
			 n = (yn - y0)/fabs(hy) + 1;
		else  n = log((yn - y0)*(ky-1)/fabs(hy) + 1)/log(ky) + 1;

		if(i == 0) flag = true;
		else flag = false;

		generate_1unreg_grid(y0,yn,hy,ky,y_grid_mas,n, flag);
	}
	sort(y_grid_mas.begin(),y_grid_mas.end());

	m_x = x_grid_mas.size();
	m_y = y_grid_mas.size();
	m_z = z_grid_mas.size();

	form_ku_for_Aker(m_x, m_y, m_z);

	printf("nodes = %d\n", m_x*m_y*m_z);
	f = fopen("xyz.txt", "w");
	fprintf(f, "%d\n", m_x*m_y*m_z);
	for(int i=0; i < m_z; i++)
		for(int j=0; j < m_y; j++)
			for(int k=0; k < m_x; k++)
				fprintf(f, "%.15lf\t%.15lf\t%.15lf\n", x_grid_mas[k], y_grid_mas[j], z_grid_mas[i]);
	fclose(f);
//кол-во ребер
	int k_x = m_x - 1;
	int k_y = m_y - 1;
	int k_z = m_z - 1 ;
	f = fopen("nvtr_3D.txt", "w");
	fprintf(f, "%d\n", (k_x )*(k_z )*(k_y));
	fprintf(f, "%d\n",(k_z + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + (k_z)*(k_x + 1)*(k_y + 1));
	printf("edges = %d\n", (k_z + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + (k_z)*(k_x + 1)*(k_y + 1));
for(int k =0; k < k_z; k++)
	for(int j = 0; j < k_y; j++)
	{
		for(int i = 0; i < k_x; i++)
		{
			//низ x
			int edge_i0 = j*k_x + j*(k_x + 1) + k*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1) + i;
			int edge_i1 = (j + 1)*k_x + (j + 1)*(k_x + 1)+ k*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1)+ i;
			//верх x
			int edge_i2 = j*k_x + j*(k_x + 1) + (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + (k+1)*(k_x + 1)*(k_y + 1) + i;
			int edge_i3 = (j + 1)*k_x + (j + 1)*(k_x + 1)+ (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + (k + 1)*(k_x + 1)*(k_y + 1)+ i;
			//низ y
			int edge_i4 = (j + 1)*k_x + j*(k_x + 1) + k*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1) + i;
			int edge_i5 = (j + 1)*k_x + j*(k_x + 1) + k*(k_x*(k_y +1) + k_y*(k_x+1)) +  k*(k_x + 1)*(k_y + 1) + i + 1;
			//верх y
			int edge_i6 = (j + 1)*k_x + j*(k_x + 1) + (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + (k + 1)*(k_x + 1)*(k_y + 1) + i;
			int edge_i7 = (j + 1)*k_x + j*(k_x + 1) + (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + (k + 1)*(k_x + 1)*(k_y + 1) + i + 1;
			//перед z
			int edge_i8 = (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
			int edge_i9 = edge_i8 + 1;
			//зад z
			int edge_i10 = (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1)+ (j + 1)*(k_x + 1) + i;
			int edge_i11 = edge_i10 + 1;

			fprintf(f, "%d %d %d %d %d %d %d %d %d %d %d %d \n", edge_i0, edge_i1, edge_i2, edge_i3, edge_i4, edge_i5, edge_i6, edge_i7, edge_i8, edge_i9, edge_i10, edge_i11);
		}
	}
	fclose(f);

	f = fopen("nvtr_v.txt", "w");
	FILE *ff = fopen("nvtr_for_tecplot.txt", "w");

	int total_elements_n = (m_z - 1) * (m_y - 1) * (m_x - 1);
	fprintf(f, "%d\n", total_elements_n);

	for(int k=0; k < m_z - 1; k++)
		for(int j=0; j < m_y - 1; j++)
			for(int i=0; i < m_x- 1; i++)
			{
				int node_i0 = k*(m_x*m_y) + j*m_x + i;
				int node_i1 = node_i0 + 1;

				int node_i2 = k*(m_x*m_y) + (j + 1)*m_x + i;
				int node_i3 = node_i2 + 1;

				int node_i4 = (k+1)*(m_x*m_y) + j*m_x + i;
				int node_i5 = node_i4 + 1;

				int node_i6 = (k+1)*(m_x*m_y) + (j + 1)*m_x + i;
				int node_i7 = node_i6 + 1;

				fprintf(ff, "%d %d %d %d %d %d %d %d\n", node_i0+1, node_i2+1, node_i4+1, node_i6+1, node_i1+1, node_i3+1, node_i5+1, node_i7+1);
				fprintf(f, "%d %d %d %d %d %d %d %d\n", node_i0, node_i1, node_i2, node_i3, node_i4, node_i5, node_i6, node_i7);

			}
	fclose(f);
	fclose(ff);
	///
	generate_ku(k_x,k_y,k_z);
	generate_ku_faces(k_x, k_y, k_z);
	return true;

	
}
bool generate_points()
{
	int n_x, n_y, n_z;//количество подотрезков
	int m_x, m_y, m_z;//количество узлов
	double x0,xn,y0,yn, z0,zn, hx, hy, hz,kx, ky,kz;
	vector<double> x_grid_mas;
	vector<double> z_grid_mas;
	vector<double> y_grid_mas;
	int n;
	bool flag = true;
	FILE *f = fopen("data_points.txt", "r");
	fscanf(f,"%d",&n_x);
	for(int i = 0; i < n_x; i++)
	{
		fscanf(f,"%lf %lf %lf %lf", &x0, &xn, &hx, &kx);
		if(x0 > xn ) return false; //проверка корректности данных
		
		if(kx == 1)
			n = (xn - x0)/fabs(hx) + 1;
		else n = log((xn - x0)*(kx-1)/fabs(hx) + 1)/log(kx) + 1;

		if(i == 0) flag = true;
		else flag = false;		

		generate_1unreg_grid(x0,xn,hx,kx,x_grid_mas,n,flag);
	}
	sort(x_grid_mas.begin(),x_grid_mas.end());

	//f = fopen("data_subseg_z.txt", "r");
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
	

	//f = fopen("data_subseg_y.txt", "r");
	fscanf(f,"%d",&n_y);
	for(int i = 0; i < n_y; i++)
	{
		fscanf(f,"%lf %lf %lf %lf", &y0, &yn, &hy, &ky);
		if( ky == 1)
			 n = (yn - y0)/fabs(hy) + 1;
		else  n = log((yn - y0)*(ky-1)/fabs(hy) + 1)/log(ky) + 1;

		if(i == 0) flag = true;
		else flag = false;

		generate_1unreg_grid(y0,yn,hy,ky,y_grid_mas,n, flag);
	}
	sort(y_grid_mas.begin(),y_grid_mas.end());

	m_x = x_grid_mas.size();
	m_y = y_grid_mas.size();
	m_z = z_grid_mas.size();

	f = fopen("points.txt", "w");
	fprintf(f, "%d\n", m_x*m_y*m_z);
	for(int i=0; i < m_z; i++)
		for(int j=0; j < m_y; j++)
			for(int k=0; k < m_x; k++)
				fprintf(f, "%lf\t%lf\t%lf\n", x_grid_mas[k], y_grid_mas[j], z_grid_mas[i]);
	fclose(f);
}
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

void matching_edges_and_nodes()
{
	FILE *nodes = fopen("nvtr_v.txt", "r");
	FILE *edges = fopen("nvtr.txt", "r");
	int N,m;


	fscanf(edges, "%d", &N);
	fscanf(edges, "%d", &m);
	KE_ed_nod *el = new KE_ed_nod[N];

	for(int i=0; i<N; i++)
	{
		for(int j=0; j<12; j++)
			fscanf(edges, "%d", &el[i].edge[j]);
	}
	fclose(edges);

		
	for(int i=0; i<N; i++)
	{
		for(int j=0; j<8; j++)
			fscanf(nodes, "%d ", &el[i].node[j]);
	}
	fclose(nodes);	

	FILE *f = fopen("nodes_edges.txt", "w");
	for(int i=0; i<N; i++)
	{
		fprintf(f, "%d %d %d\n", el[i].edge[0], el[i].node[0], el[i].node[1]);
		fprintf(f, "%d %d %d\n", el[i].edge[1], el[i].node[2], el[i].node[3]);
		fprintf(f, "%d %d %d\n", el[i].edge[2], el[i].node[4], el[i].node[5]);
		fprintf(f, "%d %d %d\n", el[i].edge[3], el[i].node[6], el[i].node[7]);
		fprintf(f, "%d %d %d\n", el[i].edge[4], el[i].node[0], el[i].node[2]);
		fprintf(f, "%d %d %d\n", el[i].edge[5], el[i].node[1], el[i].node[3]);
		fprintf(f, "%d %d %d\n", el[i].edge[6], el[i].node[4], el[i].node[6]);
		fprintf(f, "%d %d %d\n", el[i].edge[7], el[i].node[5], el[i].node[7]);
		fprintf(f, "%d %d %d\n", el[i].edge[8], el[i].node[0], el[i].node[4]);
		fprintf(f, "%d %d %d\n", el[i].edge[9], el[i].node[1], el[i].node[5]);
		fprintf(f, "%d %d %d\n", el[i].edge[10], el[i].node[2], el[i].node[6]);
		fprintf(f, "%d %d %d\n", el[i].edge[11], el[i].node[3], el[i].node[7]);
	}
	fclose(f);

	delete [] el;
}

void form_ku_for_Aker(int m_x, int m_y, int m_z)
{
	FILE *f = fopen("ku_ker.txt", "w");
	fprintf(f, "%d\n", 2*m_x*m_y + 2*(m_z - 2)*m_x + 2*(m_z - 2)*(m_y - 2));

	//нижний слой
	for(int j=0; j < m_y; j++)
		for(int i=0; i < m_x; i++)
		{
			int node_i0 = j*m_x + i;				
			fprintf(f,"%d\n", node_i0);
		}

		//верхний слой
	for(int j=0; j < m_y; j++)
		for(int i=0; i < m_x; i++)
		{
			int node_i0 = (m_z-1)*(m_x*m_y) + j*m_x + i;
			fprintf(f,"%d\n", node_i0);
		}

	//узлы по z перед
	for(int k=1; k < m_z - 1; k++)
			for(int i=0; i < m_x; i++)
			{
				int node_i4 = k*(m_x*m_y) + i;
				fprintf(f," %d\n", node_i4);
			}
	//узлы по z зад
	for(int k=1; k < m_z - 1; k++)
			for(int i=0; i < m_x; i++)
			{
				int node_i4 = k*(m_x*m_y) + (m_y - 1)*m_x + i;
				fprintf(f,"%d\n",node_i4);
			}
	//узлы по z слева
	for(int k = 1; k < m_z-1; k++)
		for(int j=1; j< m_y - 1; j++)
		{
			int node_i4 = k*(m_x*m_y) + j*m_x;
			fprintf(f,"%d\n",node_i4);
		}

	//узлы по z справа
	for(int k = 1; k < m_z-1; k++)
		for(int j=1; j< m_y-1; j++)
		{
			int node_i4 = k*(m_x*m_y) + j*m_x + m_x - 1;
			fprintf(f,"%d\n", node_i4);
		}

	fclose(f);
}

void generate_ku(int k_x, int k_y, int k_z)
{
	FILE *f = fopen("ku_3D.txt", "w");
	fprintf(f, "%d\n",2*(k_x*(k_y+1) + k_y*(k_x+1)) + 2*k_z*(k_x+1) + 2*k_z*(k_y - 1));

	//нижний слой
	for(int j=0; j<=k_y; j++)
		for(int i=0; i<k_x; i++)
		{
				int node_x = j*(k_x+1) + j*k_x + i;		
				fprintf(f, "%d %d\n", node_x, 0);
		}
	for(int j=0; j<k_y; j++)
		for(int i=0; i<=k_x; i++)
		{
				int node_y = (j+1)*(k_x) + j*(k_x + 1) + i;	
				fprintf(f, "%d %d\n", node_y,1);
		}

	//верхний слой
	for(int j=0; j<=k_y; j++)
		for(int i=0; i<k_x; i++)
		{
				int node_x = (k_z)*(k_x*(k_y +1) + k_y*(k_x+1)) + (k_z)*(k_x + 1)*(k_y + 1) + j*(k_x + 1) + j*k_x + i;		
				fprintf(f, "%d %d\n", node_x, 0);
		}
	for(int j=0; j<k_y; j++)
		for(int i=0; i<=k_x; i++)
		{
				int node_y = (k_z)*(k_x*(k_y +1) + k_y*(k_x+1)) + (k_z)*(k_x + 1)*(k_y + 1) + (j+1)*(k_x) + j*(k_x + 1) + i;	
				fprintf(f, "%d %d\n", node_y, 1);
		}

		//ребра по z перед
	for(int k=0; k < k_z; k++)
		for(int i=0; i<= k_x; i++)
		{
			int node_z = (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1) + i;
			fprintf(f, "%d %d\n", node_z, 2);
		}
		//ребра по z зад
	for(int k=0; k < k_z; k++)
		for(int i=0; i<= k_x; i++)
		{
			int node_z = (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1)+ (k_y)*(k_x + 1) + i;
			fprintf(f, "%d %d\n", node_z, 2);
		}
		//ребра по z слева
	for(int k = 0; k < k_z; k++)
		for(int j=1; j< k_y; j++)
		{
			int node_z = (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1)+  j*(k_x + 1);
			fprintf(f, "%d %d\n", node_z, 2);
		}

		//ребра по z справа
	for(int k = 0; k < k_z; k++)
		for(int j=1; j< k_y; j++)
		{
			int node_z = (k+1)*(k_x*(k_y +1) + k_y*(k_x+1)) + k*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + k_x;
			fprintf(f, "%d %d\n", node_z, 2);
		}	

	fclose(f);

	int m_x = k_x + 1;
	int m_y = k_y + 1;
	int m_z = k_z + 1;

	f = fopen("ku_v.txt", "w");
	//нижний слой
	for(int j=0; j < m_y; j++)
			for(int i=0; i < m_x - 1; i++)
			{
				int node_i0 = j*m_x + i;
				int node_i1 = node_i0 + 1;
				fprintf(f,"%d %d\n", node_i0, node_i1);
			}
	for(int j=0; j < m_y - 1; j++)
			for(int i=0; i < m_x; i++)
			{
				int node_i0 = j*m_x + i;
				int node_i2 = (j + 1)*m_x + i;
				fprintf(f,"%d %d\n", node_i0, node_i2);
			}
	//верхний слой
	for(int j=0; j < m_y; j++)
			for(int i=0; i < m_x - 1; i++)
			{
				int node_i0 = (m_z-1)*(m_x*m_y) + j*m_x + i;
				int node_i1 = node_i0 + 1;
				fprintf(f,"%d %d\n", node_i0, node_i1);
			}
	for(int j=0; j < m_y - 1; j++)
			for(int i=0; i < m_x; i++)
			{
				int node_i0 = (m_z-1)*(m_x*m_y) + j*m_x + i;
				int node_i2 = (m_z-1)*(m_x*m_y) + (j + 1)*m_x + i;
				fprintf(f,"%d %d\n", node_i0, node_i2);
			}
	//ребра по z перед
	for(int k=0; k < m_z - 1; k++)
			for(int i=0; i < m_x; i++)
			{
				int node_i0 = k*(m_x*m_y) + i;
				int node_i4 = (k + 1)*(m_x*m_y) + i;
				fprintf(f,"%d %d\n", node_i0, node_i4);
			}
	//ребра по z зад
	for(int k=0; k < m_z - 1; k++)
			for(int i=0; i < m_x; i++)
			{
				int node_i0 = k*(m_x*m_y) + (m_y - 1)*m_x + i;
				int node_i4 = (k + 1)*(m_x*m_y) + (m_y - 1)*m_x + i;
				fprintf(f,"%d %d\n", node_i0, node_i4);
			}
	//ребра по z слева
	for(int k = 0; k < m_z-1; k++)
		for(int j=1; j< m_y - 1; j++)
		{
			int node_i0 = k*(m_x*m_y) + j*m_x;
			int node_i4 = (k + 1)*(m_x*m_y) + j*m_x;
			fprintf(f,"%d %d\n", node_i0, node_i4);
		}

	//ребра по z справа
	for(int k = 0; k < m_z-1; k++)
		for(int j=1; j< m_y-1; j++)
		{
			int node_i0 = k*(m_x*m_y) + j*m_x + m_x - 1;
			int node_i4 = (k + 1)*(m_x*m_y) + j*m_x + m_x - 1;
			fprintf(f,"%d %d\n", node_i0, node_i4);
		}

	fclose(f);
}

void cut(int m_x, int m_y, int m_z, int num_kz, vector <double> &x_grid, vector <double> &y_grid, vector <double> &z_grid )//num_kz - номер слоя по z
{
	int k_x = m_x - 1;
	int k_y = m_y - 1;
	int k_z = m_z - 1;
	//выведем nvtr0
	FILE *f = fopen("nvrt0.txt", "w");
	fprintf(f, "%d\n", (m_x - 1)*(m_y - 1));
	for(int j = 0; j < m_y - 1; j++)
	{
		for(int i = 0; i < m_x - 1; i++)
		{
			int node_i0 = j*m_x + i;
			int node_i1 = j*m_x + i + 1;
			int node_i2 = (j+1)*m_x + i;
			int node_i3 = (j+1)*m_x + i + 1;

			fprintf(f, "%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
		}
	}
	
	fclose(f);

	char name[40];
	sprintf(name, "nvtr%dz.txt", num_kz);
	
	f = fopen(name, "w");
		
	for(int j=0; j < m_y - 1; j++)
	{
		for(int i=0; i < m_x - 1; i++)
		{
			int node_i0 = num_kz*(m_x*m_y) + j*m_x + i;
			int node_i1 = node_i0 + 1;

			int node_i2 = num_kz*(m_x*m_y) + (j + 1)*m_x + i;
			int node_i3 = node_i2 + 1;	
				
			fprintf(f, "%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
		}
	}
	
	fclose(f);
	int node_i0;
	int edx1, edx2, edy1, edy2, edz1, edz2;
	f = fopen("node_edges.txt", "w");
	for(int j=0; j < m_y; j++)
	{
		for(int i=0; i < m_x; i++)
		{
			node_i0 = j*m_x + i;

			if(num_kz != 0 && num_kz != m_z - 1)
			{
				if(j == 0)

					if(i == 0)
					{
						edx1 = -1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = -1;
						edy2 =  (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					}
					else if( i == m_x - 1)
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = -1;
						edy1 = -1;
						edy2 =  (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;

					}
					else
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = -1;
						edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;

					}
				else if( j == m_y - 1)
					if(i == 0)
					{
						edx1 = -1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy2 =  -1;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;

					}
					else if( i == m_x - 1)
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = -1;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)  + i ;
						edy2 =  -1;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;

					}
					else
					{
						edx1 =  j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy2 = -1;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;

					}
				else if(i == 0)
				{
					edx1 = -1;
					edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;

				}
				else if( i == m_x - 1)
				{
					edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
					edx2 = -1;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;

				}
				else
				{
					edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
					edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;

				}


			}

			else	if( num_kz == 0)
			{

				if(j == 0)

					if(i == 0)
					{
						edx1 = -1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = -1;
						edy2 =  (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = -1;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					}
					else if( i == m_x - 1)
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = -1;
						edy1 = -1;
						edy2 =  (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = -1;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;					
					}
					else
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = -1;
						edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = -1;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					}
				else if( j == m_y - 1)
					if(i == 0)
					{
						edx1 = -1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy2 =  -1;
						edz1 = -1;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					}
					else if( i == m_x - 1)
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = -1;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)  + i ;
						edy2 =  -1;
						edz1 = -1;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					}
					else
					{
						edx1 =  j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy2 = -1;
						edz1 = -1;
						edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					}
				else if(i == 0)
				{
					edx1 = -1;
					edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = -1;
					edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
				}
				else if( i == m_x - 1)
				{
					edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
					edx2 = -1;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = -1;
					edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
				}
				else
				{
					edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
					edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = -1;
					edz2 = (num_kz + 1)*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
				}
			}
			else if( num_kz == m_z - 1)
			{
				if(j == 0)

					if(i == 0)
					{
						edx1 = -1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = -1;
						edy2 =  (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = (num_kz )*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = -1;
					}
					else if( i == m_x - 1)
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = -1;
						edy1 = -1;
						edy2 =  (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = -1;

					}
					else
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = -1;
						edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edz1 = (num_kz )*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = -1;

					}
				else if( j == m_y - 1)
					if(i == 0)
					{
						edx1 = -1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy2 =  -1;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = -1;

					}
					else if( i == m_x - 1)
					{
						edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = -1;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1)  + i ;
						edy2 =  -1;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = -1;

					}
					else
					{
						edx1 =  j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
						edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
						edy2 = -1;
						edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
						edz2 = -1;

					}
				else if(i == 0)
				{
					edx1 = -1;
					edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					edz2 = -1;

				}
				else if( i == m_x - 1)
				{
					edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
					edx2 = -1;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					edz2 = -1;

				}
				else
				{
					edx1 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i - 1;
					edx2 = j*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x + 1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy1 = (j)*k_x + (j - 1)*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edy2 = (j + 1)*k_x + j*(k_x + 1) + num_kz*(k_x*(k_y +1) + k_y*(k_x+1)) + num_kz*(k_x + 1)*(k_y + 1) + i;
					edz1 = (num_kz)*(k_x*(k_y +1) + k_y*(k_x+1)) + (num_kz - 1)*(k_x + 1)*(k_y + 1)+  j*(k_x + 1) + i;
					edz2 = -1;

				}
			}
			fprintf(f, "%d\t%d %d %d %d %d %d\n", node_i0, edx1, edx2, edy1, edy2, edz1, edz2);
		}
		
	}				
		
	fclose(f);
	
	f = fopen("xyz_k.txt", "w");
	fprintf(f, "%d\n", m_x*m_y);
		for(int j=0; j < m_y; j++)
			for(int k=0; k < m_x; k++)
				fprintf(f, "%lf\t%lf\t%lf\n", x_grid[k], y_grid[j], z_grid[num_kz]);
				
	fclose(f);
	
	//system("pause");
	}

void generate_ku_faces(int k_x, int k_y, int k_z)
{

	int m_x = k_x + 1;
	int m_y = k_y + 1;
	int m_z = k_z + 1;

	FILE *f = fopen("ku_faces.txt", "w");

	int count = 2 * k_x * k_y + 2 * k_x * k_z + 2 * k_y*k_z;
	fprintf(f, "%d\n", count);

	//нижний слой
	for(int j=0; j < k_y; j++)
			for(int i=0; i < k_x; i++)
			{
				int node_i0 = j*m_x + i;
				int node_i1 = node_i0 + 1;
				int node_i2 = (j+1)*m_x + i;
				int node_i3 = node_i2 + 1;
				fprintf(f,"%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
			}
	
	//верхний слой
	for(int j=0; j < k_y; j++)
			for(int i=0; i < k_x; i++)
			{
				int node_i0 = (m_z-1)*(m_x*m_y) + j*m_x + i;
				int node_i1 = node_i0 + 1;
				int node_i2 = (m_z-1)*(m_x*m_y) + (j + 1)*m_x + i;
				int node_i3 = node_i2 + 1;
				fprintf(f,"%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
			}

	//ребра по z перед
	for(int k=0; k < k_z; k++)
			for(int i=0; i < k_x; i++)
			{
				int node_i0 = k*(m_x*m_y) + i;
				int node_i1 = node_i0 + 1;
				int node_i2 = (k + 1)*(m_x*m_y) + i;
				int node_i3 = node_i2 + 1;
				fprintf(f,"%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
			}
	//ребра по z зад
	for(int k=0; k < k_z; k++)
			for(int i=0; i < k_x; i++)
			{
				int node_i0 = k*(m_x*m_y) + (m_y - 1)*m_x + i;
				int node_i1 = node_i0 + 1;
				int node_i2 = (k + 1)*(m_x*m_y) + (m_y - 1)*m_x + i;
				int node_i3 = node_i2 + 1;
				fprintf(f,"%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
			}
	//ребра по z слева
	for(int k = 0; k < k_z; k++)
		for(int j=0; j< k_y; j++)
		{
			int node_i0 = k*(m_x*m_y) + j*m_x;
			int node_i1 = k*(m_x*m_y) + (j + 1)*m_x;
			int node_i2 = (k + 1)*(m_x*m_y) + j*m_x;
			int node_i3 = (k + 1)*(m_x*m_y) + (j + 1)*m_x;
			fprintf(f,"%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
		}

	//ребра по z справа
	for(int k = 0; k < k_z; k++)
		for(int j=0; j< k_y; j++)
		{
			int node_i0 = k*(m_x*m_y) + j*m_x + m_x - 1;
			int node_i1 = k*(m_x*m_y) + (j + 1)*m_x + m_x - 1;
			int node_i2 = (k + 1)*(m_x*m_y) + j*m_x + m_x - 1;
			int node_i3 = (k + 1)*(m_x*m_y) + (j + 1)*m_x + m_x - 1;
			fprintf(f,"%d %d %d %d\n", node_i0, node_i1, node_i2, node_i3);
		}

	fclose(f);
}
