#pragma once
#define _CRT_SECURE_NO_WARNINGS 1
#include <stdio.h>
#include <conio.h>
#include <math.h>
#include "solver_MSG.h"
#include <set>
#include "Gauss.h"

//#define NOISE



using namespace std;

struct KE //конечный элемент
{
	int area; //номер области куда попадает
	int node[4];//узлы(глобальные)
};

struct node // узел
{
	double r;//координата r
	double z;//координата z
};

struct ku_o//ребро
{
	int node; //узлы
	double ub; //значение в узлах
	
};

struct iter_data {

	double *gg;
	double *di;
	double *q0;
	double *q1;
	double *q2;
	double *q_res;
	double *pr;
	double r_con;
	double z_con;

	double *Efi0, *Efi1, *Efi2;
	double *Efi_res;

	void init();

};

class normal_field_2D
{
	private:
			
			//MSG_my s_MSG;

		void init_iter_data(iter_data* it_data);
			
			
			KE *el; //
			node *nodes;
			ku_o *ku_one;
			node *nodes_rec;
			set <int> onekuset;

			int N;//кол-во КЭ
			int m; //кол-во узлов
			int k1;//кол-во ребер с разными ку
			int n_gg;//кол-во внедиаг
			double *time_mas;
			int *ig;
			int *jg;

			int n_time;// кол-во временных слоев			
			
			double *u;
			//double *x0;
		
			void read_ku();
			void read_KE();
			void read_time();
			void read_move_controle();
			
			double form_u(double r, double z, int i, int t_i);
			double form_F(double r, double z, int i, int t_i);
			double form_mu(double r, double z, int i, int t_i);
			double form_sigma(double r, double z, int i, int t_i);
			int form_area(double r, double z1, double z2);
			void form_file_nkvat();
			
			
			void form_local_matrix_2sl(int i, int t_i, double L[4][4], double b[4], iter_data* it_data);
			void form_local_matrix_3sl(int i, int t_i, double L[4][4], double b[4], iter_data* it_data);
			void form_local_matrix_stat(int i, int t_i, double L[4][4], double b[4]);
			
			void ku(int t_i, iter_data* it_data);
			
		
			double fi_1(int i, double r, double z);
			double fi_2(int i, double r, double z);
			double fi_3(int i, double r, double z);
			double fi_4(int i, double r, double z);

			double pr_fi_1(int i, double r, double z);
			double pr_fi_2(int i, double r, double z);
			double pr_fi_3(int i, double r, double z);
			double pr_fi_4(int i, double r, double z);

			void mult_matrix(double *x, double *&y);
			void printf_matrix(int i);

			int find_KE(double r, double z);
			void Efi(int t_i, iter_data* it_data);
			double Efi_inpoint(double r, double z, iter_data* it_data);

			int mulin(int i);
			int nulin(int j);
			void form_matrix_D1(double r, double hr, double hz, double D[4][4]);
			void form_matrix_M1(double rp, double M[4][4], double hr, double hz);
			void form_matrix_G2(double hr, double hz, double r1, double r2, double G[4][4]);

			

			void solver_Gauss(double **A, double *b, double *&rez, int n);
			int n_sigma;
			double *sigma_sol, *true_sigma;
			double **rec_true, **rec1, **rec2;
			int n_rec;



			//move
			double r_con, z0_c, u_con, hcon;
			int n_con;
			double **m_rec_true;// **m_rec1;//, **m_rec2;
			
			void controle(iter_data* it_data);
			
			FILE *Efi_time;
			FILE *Efi_z;
			
			
	public:
			void read();
			void form_global_matrix(int t_i, iter_data* it_data);
			void time_calc_matrix(int i_con, double *fsigma, double **rec, iter_data* it_data);

			void calc_u(int t_i, FILE *fi, FILE *f1, double *qq);
			double Afi_inpoint(double r, double z, iter_data* it_data);
			
			void move_inverse();

			
			void move_data_true();
			int moving_controle(double *fsigma, double ** rec);
			
			
};
	