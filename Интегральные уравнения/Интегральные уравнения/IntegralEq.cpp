#include "IntegralEq.h"

const double pi = 4 * atan(1.0);
const double mu0 = 4 * pi * 1e-7;

void IntegralEq::set_env(double w_s, double sigma_0_s, double sigma_1_s, point obj_b, point obj_t) {
	w = w_s;
	sigma_0 = sigma_0_s;
	sigma_1 = sigma_1_s;
	object_bottom = obj_b;
	object_top = obj_t;

	sigma_diff = sigma_0_s - sigma_1_s;

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
	return x >= object_bottom.x && x <= object_top.x &&
		y >= object_bottom.y && y <= object_top.y &&
		z >= object_bottom.z && z <= object_top.z;

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
		nodes.push_back(node(x, y, z, (int)i));
	}
	cout << "Reading meth: reading nodes completed\n";

	dofs_n = 0;
	for(size_t el_i = 0; el_i < all_elements_n && !inp_f.eof(); el_i++) {
		vector<node> el_nodes;
		vector<dof_type> el_dofs;
		el_nodes.reserve(8);
		int tmp_int1, tmp_int2;
		inp_f >> tmp_int1 >> tmp_int2;
		for(int n_i = 0; n_i < 8; n_i++) {
			int i;
			inp_f >> i;
			i--;
			el_nodes.push_back(nodes[i]);
		}
		auto elem = new brick(el_nodes, el_dofs);
		point cent = elem->get_center();
		bool in_obj = in_object(cent.x, cent.y, cent.z);

		for (size_t n_i = 0; n_i < 8 && !in_obj; n_i++) {
			point pn = elem->get_node(n_i);
			in_obj = in_obj || in_object(pn.x, pn.y, pn.z);
		}

		if (in_obj) {
			for(int i = 0; i < 3; i++) {
				elem->add_dof(dofs_n);
				dofs_n++;
			}
			elem->set_env(w, sigma_diff);
			object_elements.push_back(elem);
		}
		all_elements.push_back(elem);
	}
	inp_f.close();
	all_elements_n = all_elements.size();
	object_elements_n = object_elements.size();

	ofstream p("obj.txt");
	for (auto& el : object_elements) {
		auto cent = el->get_center();
		p << cent.x << " " << cent.y << " " << cent.z << endl;
	}
	p.close();

}

void IntegralEq::calculate() {
	dcomplex k_val = dcomplex(1.0 / sqrt(2.0), - 1.0 / sqrt(2.0)) * sqrt(sigma_0 * w * mu0);
	G = [&, k_val](double x, double y, double z)->tfunc3d {
			return [&, x, y, z](double tx, double ty, double tz)->cmatrix(3) {
				cmatrix(3) tres;
				//vec3d dist(point(x, y, z), point(tx, ty, tz));
				vec3d dist(point(tx, ty, tz), point(x, y, z));
				double r = dist.norm();
				dcomplex exp_mult = exp(dcomplex(0, 1) * k_val * r);

				dcomplex mult1 = exp_mult / (4 * pi * r * r * r * sigma_0);
				dcomplex common_mult = mult1 * (-k_val * k_val * r * r + 3.0 * dcomplex(0, 1) * k_val * r + 3.0);
				dcomplex k_sq = k_val * k_val;
				
				for(int i = 0; i < 3; i++) {
					for (int j = 0; j < 3; j++) {
						if (i == j)
							tres[i][j] = mult1 * ((-k_sq * r * r + 3.0 * dcomplex(0, 1) * k_val * r + 3.0) * dist[i] * dist[j] / (r * r) + (k_sq * r * r - dcomplex(0, 1) * k_val * r - 1.0));
						else
							tres[i][j] = mult1 * (-k_sq * r * r + 3.0 * dcomplex(0, 1) * k_val * r + 3.0) * dist[i] * dist[j] / (r * r);
					}
					//tres[i][i] = tres[i][i] + mult1 * (-k_val * k_val * r * r - dcomplex(0, 1) * k_val * r - 1.0);
 				}

				return tres;
			};
		
	};

	dyn_matrix A_glob;
	A_glob.resize(dofs_n);
	for(size_t i = 0; i < dofs_n; i++)
		A_glob[i].resize(dofs_n, 0);
	vector<dcomplex> b_glob;
	b_glob.resize(dofs_n, 0);

	cout << "Calculate: assambeling matrix\n";

	for(size_t el_i = 0; el_i < object_elements_n; el_i++) {

		point el_cent = object_elements[el_i]->get_center();
		tfunc3d G_p = G(el_cent.x, el_cent.y, el_cent.z);

		auto row_dofs = object_elements[el_i]->get_dofs_num();

		for(auto&el_E0 : object_elements) {

			// Расчёт вклада в матрицу
			cmatrix(3) el_matrix = el_E0->get_matrix_value(G_p, el_cent);
			auto col_dofs = el_E0->get_dofs_num();
			for(size_t a_i = 0; a_i < row_dofs.size(); a_i++) {
				for(size_t a_j = 0; a_j < col_dofs.size(); a_j++)
					A_glob[row_dofs[a_i]][col_dofs[a_j]] += el_matrix[a_i][a_j];
			}


			// Расчёт вклада в правую часть
			vfunc3d el_integ_func = [&](double x, double y, double z)->vec3d {
				vec3d E0_val = calc_E0(x, y, z);
				cmatrix(3) G_val = G_p(x, y, z);
				vec3d integ_res = w * (sigma_diff) * (G_val * E0_val);
				return integ_res;
			};
			vec3d el_add = el_E0->integrate(el_integ_func);

			for(size_t b_i = 0; b_i < row_dofs.size(); b_i++)
				b_glob[row_dofs[b_i]] += el_add[b_i];
		}

	}


	Ep_solution.clear();
	Ep_solution.resize(dofs_n);

	ofstream tt("tt.txt");
	for (int i = 0; i < 3 * object_elements_n; i++) {
		for (int j = 0; j < 3 * object_elements_n; j++) {
			tt << A_glob[i][j] << "\t";
		}
		tt << endl;
	}
	tt.close();

	cout << "Calculate: solving\n";
	SLAE_solution_Gauss(A_glob, b_glob.data(), Ep_solution.data(), (int)dofs_n);
	cout << "Calculate: solving completed\n";
}

vec3d IntegralEq::calc_E0(double x, double y, double z) {
	vec3d res(0, 0, 0);

	tfunc3d G_p = G(x, y, z);

	for(auto& sect : J_sect) {
		vec3d sect_res = sect->integ_dir(G_p);
		res = res + sect_res;
	}

	double v1 = 0;
	for (auto& sect : J_sect) {
		func3d l = [](double x, double y, double z) {return 1.0 / 4.0; };
		v1 += sect->integrate(l);
	}

	return res;
}

vec3d IntegralEq::calc_Ep(double x, double y, double z) {
	vec3d res(0, 0, 0);

	if (in_object(x, y, z)) {
		bool found = false;
		for (auto& el : object_elements) {
			if (el->in_element(x, y, z)) {
				auto loc_dofs = el->get_dofs();
				for (size_t dof_i = 0; dof_i < loc_dofs.size(); dof_i++) {
					res = res + Ep_solution[loc_dofs[dof_i]] * el->basis(dof_i, x, y, z);
				}
				return res;
			}
		}
	}
	tfunc3d G_p = G(x, y, z);
	for (auto& el : object_elements) {
		auto el_dofs = el->get_dofs_num();
		vfunc3d integ_f = [&](double tx, double ty, double tz)->vec3d {
			cmatrix(3) G_val = G_p(tx, ty, tz);
			vec3d E_val(0, 0, 0);
			for (int i = 0; i < el_dofs.size(); i++) {
				E_val = E_val + Ep_solution[el_dofs[i]] * el->basis(i, tx, ty, tz);
			}
			vec3d E0_val = calc_E0(x, y, z);
			E_val = E_val + E0_val;

			vec3d res = G_val * E_val;
			res = dcomplex(0, 1) * w * sigma_diff * res;
			return res;
		};
		vec3d el_res = el->integrate(integ_f);

		res = res + el_res;

	}

	return res;

}

vec3d IntegralEq::calc_E(double x, double y, double z) {
	vec3d res = calc_E0(x, y, z);

	res = res + calc_Ep(x, y, z);

	return res;

}

void IntegralEq::output(string file_name) {
	ofstream outpf(file_name.c_str());

	outpf << "VARIABLES = \"x\" \"y\" \"z\" \"AxR\" \"AyR\" \"AzR\" \"AxI\" \"AyI\" \"AzI\"\n";  

	for(auto& el : all_elements) {
		point cent = el->get_center();
		vec3d val = calc_E(cent.x, cent.y, cent.z);

		outpf << cent.x << " " << cent.y << " " << cent.z << " " << val.x.real() << " " << val.y.real() << " " << val.z.real()<< " " << val.x.imag() << " " << val.y.imag() << " " << val.z.imag() << endl;

	}
	outpf.close();
}

void IntegralEq::outputE0(string file_name) {
	ofstream outpf(file_name.c_str());

	outpf << "VARIABLES = \"x\" \"y\" \"z\" \"AxR\" \"AyR\" \"AzR\" \"AxI\" \"AyI\" \"AzI\"\n";  

	for(auto& el : all_elements) {
		point cent = el->get_center();
		vec3d val = calc_E0(cent.x, cent.y, cent.z);

		outpf << cent.x << " " << cent.y << " " << cent.z << " " << val.x.real() << " " << val.y.real() << " " << val.z.real()<< " " << val.x.imag() << " " << val.y.imag() << " " << val.z.imag() << endl;

	}
	outpf.close();
}

void IntegralEq::outputE_surface(string file_name) {
	ofstream outpf(file_name.c_str());

	outpf << "VARIABLES = \"x\" \"y\" \"z\" \"AxR\" \"AyR\" \"AzR\" \"AxI\" \"AyI\" \"AzI\"\n";  

	size_t n_x = 50, n_y = 50;
	double x_min = -50, x_max = 0, y_min = -50, y_max = 0;
	double z = 0;
	double hx = (x_max - x_min) / n_x;
	double hy = (y_max - y_min) / n_y;

	for(size_t i_x = 0; i_x < n_x; i_x++) {
		for(size_t i_y = 0; i_y < n_y; i_y++) {
			point cent(x_min + i_x * hx, y_min + i_y * hy, z);
			vec3d val = calc_E(cent.x, cent.y, cent.z);

			outpf << cent.x << " " << cent.y << " " << cent.z << " " << val.x.real() << " " << val.y.real() << " " << val.z.real()<< " " << val.x.imag() << " " << val.y.imag() << " " << val.z.imag() << endl;

		}

	}
	

	outpf.close();
}


void IntegralEq::outputE_surface_tp(string file_name) {
	ofstream outpf(file_name.c_str());

	outpf << "TITLE = \"Slice\" \n";
	outpf << "VARIABLES = \"x\" \"y\" \"AxR\" \"AyR\" \"AzR\" \"AxI\" \"AyI\" \"AzI\" \"Ax\" \"Ay\"\n";

	size_t n_x = 50, n_y = 50;
	double x_min = -5, x_max = 5, y_min = -5, y_max = 5;
	double z = -3;
	double hx = (x_max - x_min) / n_x;
	double hy = (y_max - y_min) / n_y;

	outpf << "ZONE I= " << n_x << ", J= " << n_y << ", F=POINT\n";

	double t = 440;


	for (size_t i_x = 0; i_x < n_x; i_x++) {
		for (size_t i_y = 0; i_y < n_y; i_y++) {
			point cent(x_min + i_x * hx, y_min + i_y * hy, z);
			vec3d val = calc_E(cent.x, cent.y, cent.z);

			vec3d A_real;
			dcomplex exp_val = exp(dcomplex(0, 1) * w * t);
			vec3d real_tmp = exp_val * val;
			for (int a_i = 0; a_i < 3; a_i++)
				A_real[a_i] = real_tmp[a_i].real() + real_tmp[a_i].imag();
			//for (int a_i = 0; a_i < 3; a_i++)
			//	A_real[a_i] = (val[a_i].real() * cos(w * t) - val[a_i].imag() * sin(w * t)) + (val[a_i].real() * sin(w * t) + val[a_i].imag() * cos(w * t));

			outpf << cent.x << " " << cent.y << " " << val.x.real() << " " << val.y.real() << " " << val.z.real() << " " << val.x.imag() << " " << val.y.imag() << " " << val.z.imag()  << " " << A_real.x.real() << " " << A_real.y.real() << endl;

		}

	}
	outpf << endl;

	outpf.close();
}
