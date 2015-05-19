#pragma once

#include <stdio.h>
#include <conio.h>
#include "math.h"
#include <vector>
#include <algorithm>
using namespace std;

struct KE_ed_nod //конечный элемент
{
	int node[8];//узлы(глобальные)
	int edge[12];//ребра глобальные
};
void generate_1reg_grid(double a, double b, double h, double *grid_mass, int n);
void generate_1unreg_grid(double a, double b, double h, double kr, vector<double>&grid_mass, int n, bool flag);
bool generate_2reg_grid();
bool generate_2unreg_grid();
bool generate_3unreg_grid();
void matching_edges_and_nodes();
bool generate_points();
void form_ku_for_Aker(int m_x, int m_y, int m_z);
void cut(int m_x, int m_y,int m_z, int num_kz, vector <double> &x_grid, vector <double> &y_grid, vector <double> &z_grid );//num_kz - номер слоя по z
void generate_ku(int k_x, int k_y, int k_z);
void generate_ku_faces(int k_x, int k_y, int k_z);
