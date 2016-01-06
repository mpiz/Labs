#include "elements_classes.h"
#include "Gauss.h"
#include <iostream>

class IntegralEq {
public:
	void set_env(double w_s, double sigma_0_s, double sigma_1_s, point obj_b, point obj_t);
	void input_mesh(string file_name);
	
	vec3d calc_E0(double x, double y, double z);
	vec3d calc_Ep(double x, double y, double z);
	vec3d calc_E(double x, double y, double z);
	void calculate();

	void output(string file_name);
	void outputE0(string file_name);
	void outputE_surface(string file_name);
	void outputE_surface_tp(string file_name);

private:

	size_t all_elements_n;
	size_t object_elements_n;
	size_t nodes_n;

	vector<brick*> all_elements;
	vector<brick*> object_elements;
	vector<node> nodes;
	
	double w, sigma_0, sigma_1;
	double sigma_diff;
	point object_bottom, object_top;

	bool in_object(double x, double y, double z);

	plane J_plane;

	vector<sector*> J_sect;

	function<tfunc3d(double x, double y, double z)> G;

	vector<dcomplex> Ep_solution;
	size_t dofs_n;

};