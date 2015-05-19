#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "geom_classes.h"
#include "integration_points.h"
#include <functional>
#include <cstdio>
#include <vector>
#include <algorithm>

typedef function<double(double, double, double)> func3d;
typedef function<vec3d(double, double, double)> vfunc3d;

const int gauss_points_tr = 3;
const int gauss_points_tet = 4;

class quadelement;


class quadelement {
 public:
	 quadelement();
	 quadelement(vector<node> nodes_s, vector<dof_type> s_dofs);

	  dof_type& operator [] (int i);
	  double integrate(func3d integ_func);

	  dyn_matrix get_local_matrix(double mu);
	  vector<double> get_local_right_part(vfunc3d rp_func);

	  vec3d vector_basis_v(int i, double x, double y, double z);

	  vector<dof_type> get_dofs();

	  void renumerate_dofs(map<dof_type, dof_type>& glob_to_loc);

 private:

	 void init_coords();

	 static const int element_nodes = 4;

	 unsigned int dofs_number;
	 array<node, element_nodes> node_array;
	 array<point, element_nodes> node_array_local;

	 vector<dof_type> dofs;


	 vec3d normal_vector;
	 array<vec3d, 3> tau;
	 point to_local_cord(point p_glob);
	 point to_global_cord(point p_loc);

	 matrix(3) transition_matrix; //������� �������� � ��������� ����������

	  point transform(point pr); //��������� ���������� ���������� � ���������
	 
	 double h_u, h_v;
	 double u_0, v_0;

	 double get_u(double x, double y, double z);
	 double get_v(double x, double y, double z);

	 vec3d mull_transT(vec3d v);

	 int gauss_points_n;
	 vector<point> gauss_points; //����� ��� �������������� �� ������
	 vector<point> gauss_points_global;
	 vector<double> gauss_weights;
	 double jacobian; //������� ��� ��������� ���������

};


class cubeelement {
 public:
	 cubeelement();
	 cubeelement(vector<node> nodes_s, vector<dof_type> s_dofs);

	  dof_type& operator [] (int i);
	   double integrate(func3d integ_func);

	  dyn_matrix get_local_matrix(double mu);
	  vector<double> get_local_right_part(vfunc3d rp_func);

	  vec3d vector_basis_v(int i, double x, double y, double z);
	  vec3d vector_basis_rot_v(int i, double x, double y, double z);

	  vector<dof_type> get_dofs();

	  bool in_element(point pn);
 private:

	 void init_coords();
	 static const int element_nodes = 8;

	 unsigned int dofs_number;
	 array<node, element_nodes> node_array;

	 vector<dof_type> dofs;



	 point to_local_cord(point p_glob);
	 point to_global_cord(point p_loc);

	  point transform(point pr); //��������� ���������� ���������� � ���������
	 
	 double h_x, h_y, h_z;
	 double x_0, y_0, z_0;

	 int gauss_points_n;
	 vector<point> gauss_points; //����� ��� �������������� �� ������
	 vector<point> gauss_points_global;
	 vector<double> gauss_weights;
	 double jacobian; //������� ��� ��������� ���������

};
