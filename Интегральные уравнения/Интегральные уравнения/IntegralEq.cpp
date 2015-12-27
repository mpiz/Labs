#include "IntegralEq.h"

const double pi = 4 * atan(1.0);

void IntegralEq::set_env(double w_s, double sigma_0_s, double sigma_1_s, point obj_b, point obj_t) {
	w = w_s;
	sigma_0 = sigma_0_s;
	sigma_1 = sigma_1_s;
	object_bottom = obj_b;
	object_top = object_top;

	array<node, 3> J_array_node = {{node(0, 0, 0), node(1, 0, 0), node(0, 1, 0)}};

	J_plane = plane(J_array_node);
	vector<node> J_nodes;
	J_nodes.resize(2);
	J_nodes[0] = node(-1, -1, 0, 0); J_nodes[1] = node(1, -1, 0, 1);
	J_sect.push_back(new sector(J_nodes, J_plane));

	J_nodes[0] = node(1, -1, 0, 0); J_nodes[1] = node(1, 1, 0, 1);
	J_sect.push_back(new sector(J_nodes, J_plane));

	J_nodes[0] = node(1, 1, 0, 0); J_nodes[1] = node(-1, 1, 0, 1);
	J_sect.push_back(new sector(J_nodes, J_plane));

	J_nodes[0] = node(-1, 1, 0, 0); J_nodes[1] = node(-1, -1, 0, 1);
	J_sect.push_back(new sector(J_nodes, J_plane));

}

bool IntegralEq::in_object(double x, double y, double z) {
	return x >= -100 && x <= 100 &&
		y >= -100 && y <= 100 &&
		z >= -10 && z <= -2;

}

void IntegralEq::input_mesh(string file_name) {

	ifstream inp_f(file_name.c_str());

	inp_f >> nodes_n >> all_elements_n;
	nodes.reserve(nodes_n);
	for(size_t n_i = 0; n_i < nodes_n; n_i++) {
		size_t i;
		double x, y, z;
		inp_f >> i >> x >> y >> z;
		i--;
		nodes.push_back(node(x, y, z, i));
	}

	dofs_n = 0;
	for(size_t el_i = 0; el_i < all_elements_n; el_i++) {
		vector<node> el_nodes;
		vector<dof_type> el_dofs;
		el_nodes.reserve(8);
		for(int n_i = 0; n_i < 8; n_i++) {
			int i;
			inp_f >> i;
			i--;
			el_nodes.push_back(nodes[i]);
		}
		auto elem = new brick(el_nodes, el_dofs);
		point cent = elem->get_center();
		if (in_object(cent.x, cent.y, cent.z)) {
			for(int i = 0; i < 3; i++) {
				elem->add_dof(dofs_n);
				dofs_n++;
			}
			object_elements.push_back(elem);
		}
		all_elements.push_back(elem);
	}

	object_elements_n = object_elements.size();

	G = [&](double x, double y, double z)->tfunc3d {
		dcomplex k_val = (1.0 / sqrt(2.0), - 1.0 / sqrt(2.0)) * sqrt(sigma_0 * w);
			return [&](double tx, double ty, double tz)->cmatrix(3) {
				cmatrix(3) tres;
				vec3d dist(point(x, y, z), point(tx, ty, tz));
				double r = dist.norm();
				dcomplex exp_mult = (cos(k_val * r), sin(k_val * r));

				dcomplex common_mult = exp_mult * (-k_val * k_val * r * r + 3.0 * dcomplex(0, 1) * k_val*r) / (4 * pi * r * r * r);
				
				for(int i = 0; i < 3; i++)
					for(int j = 0; j < 3; j++)
						tres[i][j] = common_mult * dist[i] * dist[j];

				return tres;
			};
		
	};
}

void IntegralEq::calculate() {

	dyn_matrix A_glob;
	A_glob.resize(dofs_n);
	for(size_t i = 0; i < dofs_n; i++)
		A_glob[i].resize(dofs_n, 0);
	vector<dcomplex> b_glob;
	b_glob.resize(dofs_n, 0);


	for(size_t el_i = 0; el_i < object_elements_n; el_i++) {
		point el_cent = object_elements[el_i]->get_center();
		tfunc3d G_p = G(el_cent.x, el_cent.y, el_cent.z);
		cmatrix(3) el_matrix = object_elements[el_i]->get_matrix_value(G_p, el_cent);

		vec3d b_loc(0, 0, 0);
		for(auto&el_E0 : object_elements) {
			vfunc3d el_integ_func = [&](double x, double y, double z)->vec3d {
				vec3d E0_val = calc_E0(x, y, z);
				vec3d integ_res = G_p(x, y, z) * E0_val;
				return integ_res;
			};
			vec3d el_add = el_E0->integrate(el_integ_func);
			b_loc = b_loc + el_add;
		}

		size_t add_i = 3 * el_i;
		for(size_t a_i = 0; a_i < 3; a_i++) {
			for(size_t a_j = 0; a_j < 3; a_j++) {
				A_glob[add_i + a_i][add_i + a_j] += el_matrix[a_i][a_j];
			}
			b_glob[add_i + a_i] += b_loc[a_i];
		}

	}
	SLAE_solution_Gauss(A_glob, b_glob.data(), Ep_solution.data(), dofs_n);
}

vec3d IntegralEq::calc_E0(double x, double y, double z) {
	vec3d res(0, 0, 0);

	tfunc3d G_p = G(x, y, z);

	for(auto& sect : J_sect) {
		vec3d sect_res = sect->integ_dir(G_p);
		res = res + sect_res;
	}

	return res;
}