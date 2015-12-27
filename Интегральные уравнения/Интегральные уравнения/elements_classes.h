#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include "add_geom_classes.h"
#include <functional>
#include <cstdio>
#include <vector>

typedef function<double(double, double, double)> func3d;
typedef function<vec3d(double, double, double)> vfunc3d;
typedef function<cmatrix(3)(double, double, double)> tfunc3d;

const int gauss_points_sec = 3;
const int gauss_points_tr = 3;
const int gauss_points_tet = 4;

class simple_element;



class simple_element {

public:
	virtual vfunc3d get_vector_basis_dof(size_t dof_i);
	virtual vfunc3d get_vector_basis(dof_type order, dof_type num = 0);


	virtual vfunc3d get_vector_right_part_dof(size_t dof_i);
	virtual vfunc3d get_vector_right_part(dof_type order, dof_type num = 0);

	virtual dyn_matrix get_local_matrix(double mu);
	virtual vector<double> get_local_right_part(func3d rp_func);
	virtual vector<double> get_local_right_part(vfunc3d rp_func);
	virtual double integrate(func3d func);	// Вычисление интеграла по элементу
	virtual vec3d integrate(vfunc3d func);	// Вычисление интеграла по элементу


	void add_dof(dof_info d);
	void prepare_gauss(int gn);
	vector<dof_info> get_dofs();
	vector<dof_type> get_dofs_num();


protected:
	
	vector<dof_info> dofs;
	unsigned int dofs_number;


	size_t gauss_points_n;
	double jacobian;
	vector<point> gauss_points_global; //точки для интегрирования по Гауссу (в локальной системе координат)
	vector<double> gauss_weights;
};


class sector : public simple_element {
public:


	sector();
	sector(const vector<node>& nodes_s, const plane& plane_s);
	sector(vector<node> nodes_s, vector<dof_info> s_dofs);

	int& operator [] (int i); //получить i-ую степень свободы

	static const int element_nodes = 2;

	bool in_element(double x, double y, double z);

	void for_point_on_element(function<void(double, double, double)> func);

	vec3d integ_dir(tfunc3d G);

private:
	vector<node> nodes;
	plane sector_plane;	// Плоскость, в которой лежит отрезок

	vec3d direction;
	double length;

	vec3d normal_in_plane;	// Нормальный вектор в плоскости sector_plane

	double get_t(double x, double y, double z);
	point get_point(double t);

	void init_coords();

};

class brick : public simple_element {
public: 
	brick(vector<node>& nodes_s, vector<dof_type> dofs_s); 

	void set_env(double w_s, double sigma_s);

	cmatrix(3) get_matrix_value(tfunc3d G, point pn);
	vec3d basis(dof_type loc_dof, double x, double y, double z);

	bool in_element(double x, double y, double z);
	point get_center();
private:
	vector<node> nodes;

	double hx, hy, hz;


	double w, sigma;

};