#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "geom_classes.h"
#include "L_cords_gen.h"
#include <functional>
#include <cstdio>

typedef function<dcomplex(double, double, double)> func3d;
typedef function<vec3d(double, double, double)> vfunc3d;

const int gauss_points_tr = 3;
const int gauss_points_tet = 4;

class trelement;
class tetelement;

typedef vec3d (trelement::*tr_bfunc)(point);
typedef vec3d (tetelement::*tet_bfunc)(point);

class trelement {
 public:
	 trelement();
	 trelement(node n1, node n2, node n3);

	 int& operator [] (int i); //�������� i-�� ��������� �����

	 double integral(func3d func); //���������� ��������� �� ������������ (�����, 3 �����) � ��������� �����������

	 int get_ph_area();
	 void set_ph_area(int sph_area);

	 matrix(6) get_M_matrix();
	 matrix(6) get_kernel_M_matrix();

	 array<dcomplex, 6> get_right_part(vfunc3d rp_func);
	 array<dcomplex, 6> get_kernel_right_part(vfunc3d rp_func);

	 array<int, 3> loc_edge;

	 dcomplex integrate(func3d integ_func);//���������� ��������� �� ������������ � ���������� �����������

	 vec3d basis_v(int i, double x, double y, double z);
	 node local_node(int i);
	 vec3d get_tau(int i);


	 
 private:

	 void generate_L_cords();

	 point to_local_cord(point p_glob);
	 point to_global_cord(point p_loc);

	 vec3d to_global_cord(vec3d v_loc);

	 double lambda(int l_i, point p_loc);

	 tr_bfunc basis[6];
	 tr_bfunc kernel_basis[6];

	 matrix(3) D_matrix, L_cord_matrix; //L - ������� L-���������, D - � ���������
	 double det_D; //������������� ������� L-���������

	 void init_cords(); //���������� ��������� ���������

	  array<node, 3> node_array;
	  array<int, 3> edge_array;

	  vec3d w1(point p_loc);
	  vec3d w2(point p_loc);
	  vec3d w3(point p_loc);
	  vec3d w4(point p_loc);
	  vec3d w5(point p_loc);
	  vec3d w6(point p_loc);

	 /* ��������� �����
	   w1 = l1 * grad(l2) - l2 * grad(l1);
	  w2 = l1 * grad(l3) - l3 * grad(l2);
	  w3 = l2 * grad(l3) - l3 * grad(l2);

	  w4 = l1 * grad(l2) + l2 * grad(l1);
	  w5 = l1 * grad(l3) + l3 * grad(l2);
	  w6 = l2 * grad(l3) + l3 * grad(l2);
	  */

	  vec3d gradphi1(point p_loc);
	  vec3d gradphi2(point p_loc);
	  vec3d gradphi3(point p_loc);

	  vec3d grad_lambda(int i);

	  matrix(3) transition_matrix; //������� �������� � ��������� ����������

	  point transform(point pr); //��������� ���������� ���������� � ���������

	  int ph_area; //���������� �������
	  point trpoint[3]; //��������� ���������� ����� ������������
	  point gauss_points[gauss_points_tr]; //����� ��� �������������� �� ������ (� ��������� ������� ���������)
	  point gauss_points_global[gauss_points_tr]; //����� ��� �������������� �� ������ (� ��������� ������� ���������)
	  double gauss_weights[gauss_points_tr];
	  double jacobian; //������� ��� ��������� ���������

	  array<vec3d, 4> tau;
};


class tetelement {
 public:
	 tetelement();
	 tetelement(node n1, node n2, node n3, node n4);

	 int& operator [] (int i); //�������� i-�� ��������� �����

	 cmatrix(12) get_local_matrix(double mu, dcomplex k_sq);
	 cmatrix(10) get_local_kernel_matrix(dcomplex k_sq);
	 matrix(12) get_local_M();
	 recmatrix(10,12) get_local_R();
	 array<dcomplex, 12> get_local_right_part(vfunc3d rp_func);

	 dcomplex integrate(func3d integ_func);

	 int get_ph_area();
	 void set_ph_area(int sph_area);

	 node get_local_node(int i);
	 point get_center();
	 vec3d basis_v(int i, double x, double y, double z);
	 vec3d basis_rot_v(int i, double x, double y, double z);

	 bool in_element(double x, double y, double z);
	 bool valid_for_tree_node(double x0, double x1, double y0, double y1, double z0, double z1);


 private:

	 tet_bfunc basis[12];
	 tet_bfunc basis_rotors[6];
	 tet_bfunc kernel_basis[10];

	 void init_coords();
	 void generate_L_cords();
	 void calculate_M_matrix();

	 double lambda(int i, point p_glob);	//�������� i-� L-���������� � �����
	 vec3d grad_lambda(int i);	//�������� i-� L-����������

	 //�������� �������� �������
	 vec3d w1(point p_glob);	//0, 1
	 vec3d w2(point p_glob);	//0, 2
	 vec3d w3(point p_glob);	//0, 3
	 vec3d w4(point p_glob);	//1, 2
	 vec3d w5(point p_glob);	//1, 3
	 vec3d w6(point p_glob);	//2, 3

	 //���������� �������� �������
	 vec3d w7(point p_glob);
	 vec3d w8(point p_glob);
	 vec3d w9(point p_glob);
	 vec3d w10(point p_glob);
	 vec3d w11(point p_glob);
	 vec3d w12(point p_glob);

	 //������ �������
	 vec3d rotw1(point p_glob);
	 vec3d rotw2(point p_glob);
	 vec3d rotw3(point p_glob);
	 vec3d rotw4(point p_glob);
	 vec3d rotw5(point p_glob);
	 vec3d rotw6(point p_glob);

	 //��������� ������� ������� ������� �������
	 vec3d gradphi1(point p_glob);
	 vec3d gradphi2(point p_glob);
	 vec3d gradphi3(point p_glob);
	 vec3d gradphi4(point p_glob);

	 array<node,4> node_array;
	 array<int,6> edge_array;

	 int ph_area; //���������� �������


	 matrix(4) D_matrix, L_cord_matrix; //L - ������� L-���������, D - � ���������
	 //i - ������ �� i-� ����������, ���� ������������ x,y,z,1
	 double det_D; //������������� ������� L-���������
	 matrix(12) M_matrix;

	  point gauss_points[gauss_points_tet]; //����� ��� �������������� �� ������
	  double gauss_weights[gauss_points_tet];
	  double jacobian; //������� ��� ��������� ���������

	  //��� ���������� ������
	  array<double, 3> ch_points[5];
	  double edges_a[6][3], edges_b[6][3];
};