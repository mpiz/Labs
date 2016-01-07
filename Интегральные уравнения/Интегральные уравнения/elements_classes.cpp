#include "elements_classes.h"

vector<dof_info> simple_element::get_dofs() {
	return dofs;
}

vector<dof_type> simple_element::get_dofs_num() {
	vector<dof_type> res;
	for(auto& dof_it : dofs) {
		res.push_back(dof_it);
	}
	return res;
}

vfunc3d simple_element::get_vector_basis_dof(size_t dof_i) {
	return get_vector_basis(1, 0);
}

vfunc3d simple_element::get_vector_basis(dof_type order, dof_type num) {
	return nullptr;
}

void simple_element::add_dof(dof_info d) {
	dofs.push_back(d);
	dofs_number = dofs.size();
}

vfunc3d simple_element::get_vector_right_part_dof(size_t dof_i) {
	return get_vector_right_part(1, 0);
}

vfunc3d simple_element::get_vector_right_part(dof_type order, dof_type num) {
	return nullptr;
}

vector<double> simple_element::get_local_right_part(func3d rp_func) {
	throw;
}

vector<double> simple_element::get_local_right_part(vfunc3d rp_func) {
	throw;
}

dyn_matrix simple_element::get_local_matrix(double mu) {
	throw;
}


double simple_element::integrate(func3d func) {

	double res = 0;

	for(size_t i = 0; i < gauss_points_n; i++) {
		res += gauss_weights[i] * func(gauss_points_global[i].x, gauss_points_global[i].y, gauss_points_global[i].z);
	}

	res *= jacobian;

	return res;
}

vec3d simple_element::integrate(vfunc3d func) {

	vec3d res(0, 0, 0);

	for(size_t i = 0; i < gauss_points_n; i++) {
		res = res + gauss_weights[i] * func(gauss_points_global[i].x, gauss_points_global[i].y, gauss_points_global[i].z);
	}

	res = jacobian * res;

	return res;
}

void simple_element::prepare_gauss(int gn) {
	gauss_points_n = gn;
	gauss_points_global.resize(gauss_points_n);
	gauss_weights.resize(gauss_points_n);

}

// ========  Отрезки ========
sector::sector() {
	prepare_gauss(gauss_points_sec);
}

sector::sector(const vector<node>& nodes_s, const plane& plane_s) {
	nodes = nodes_s;
	sector_plane = plane_s;

	// Зададим строгий порядок узлов
	if (nodes[0].number > nodes[1].number)
		swap(nodes[0], nodes[1]);

	init_coords();
}

sector::sector(vector<node> nodes_s, vector<dof_info> s_dofs) {
	dofs = s_dofs;
	dofs_number = dofs.size();
	nodes = nodes_s;

	// Зададим строгий порядок узлов
	if (nodes[0].number > nodes[1].number)
		swap(nodes[0], nodes[1]);

	init_coords();

}

void sector::init_coords() {
	if (nodes.size() != 2) {
		throw "More or less then 2 edge nodes!";
	}

	prepare_gauss(gauss_points_sec);

	direction = vec3d(nodes[0], nodes[1]);
	length = direction.norm();
	direction = direction / length;

	jacobian = length / 2.0;

	double gauss_coeff = sqrt(3.0 / 5.0);

	gauss_points_global[0] = get_point((-gauss_coeff + 1) / 2.0);
	gauss_points_global[1] = get_point(0.5);
	gauss_points_global[2] = get_point((-gauss_coeff + 1) / 2.0);

	gauss_weights[0] = 5.0 / 9.0;
	gauss_weights[1] = 8.0 / 9.0;
	gauss_weights[2] = 5.0 / 9.0;

	if (sector_plane.get_jacobian() != 0) {
		// Построим нормаль
		normal_in_plane = sector_plane.get_normal_in_plane(direction);
		normal_in_plane = normal_in_plane / normal_in_plane.norm();
	}
	else {
		normal_in_plane = vec3d(0, 0, 0);
	}
}

point sector::get_point(double t) {
	point res = (vec3d(nodes[0]) +  t * length * direction).to_point();

	return res;

}

vec3d sector::integ_dir(tfunc3d G) {
	vfunc3d G_dir = [&](double x, double y, double z)->vec3d {
		cmatrix(3) G_val = G(x, y, z);
		vec3d res = G_val * direction;
		return res;
	};

	vec3d integ_val = integrate(G_dir);

	return integ_val;
}


brick::brick(vector<node>& nodes_s, vector<dof_type> dofs_s) {
	nodes  = nodes_s;
	dofs = dofs_s;


	// Упорядочим вершины, в нужном виде
	sort(nodes.begin(), nodes.end(), [](const node& n1, const node& n2)->bool {
		return n1.z < n2.z;
	});
	sort(nodes.begin(), nodes.begin() + 4, [](const node& n1, const node& n2)->bool {
		return n1.y < n2.y;
	});

	sort(nodes.begin()+4, nodes.end(), [](const node& n1, const node& n2)->bool {
		return n1.y < n2.y;
	});

	for(auto i = 0; i < nodes.size(); i+=2) {
		if (nodes[i].x < nodes[i+1].x)
			swap(nodes[i], nodes[i+1]);
	}

	// Подготовим данные, для интегрирования
	hx = nodes[7].x - nodes[0].x;
	hy = nodes[7].y - nodes[0].y;
	hz = nodes[7].z - nodes[0].z;

	jacobian = hx * hy * hz / 8.0;
	prepare_gauss(8);

	auto trans_from_master = [&](double ksi, double etta, double dzeta)->point {
		point pn (
			hx * (ksi + 1) / 2.0 + nodes[0].x,
			hy * (etta + 1) / 2.0 + nodes[0].y,
			hz * (dzeta + 1) / 2.0 + nodes[0].z
		);
		return pn;
	};

	prepare_gauss(gauss_points_brick);

	double g_a = sqrt(6.0 / 7.0);
	double g_b = sqrt((960.0 - 33.0 * sqrt(238.0)) / 2726.0);
	double g_c = sqrt((960.0 - 33.0 * sqrt(238.0)) / 2726.0);
	double w_1 = 1078.0 / 3645.0;
	double w_2 = 343.0 / 3645.0;
	double w_3 = 43.0 / 135.0 + 829.0 * sqrt(238.0) / 136323.0;
	double w_4 = 43.0 / 135.0 - 829.0 * sqrt(238.0) / 136323.0;

	gauss_points_global[0] = trans_from_master(g_a, 0, 0);
	gauss_points_global[1] = trans_from_master(-g_a, 0, 0);

	gauss_points_global[2] = trans_from_master(0, g_a, 0);
	gauss_points_global[3] = trans_from_master(0, -g_a, 0);

	gauss_points_global[4] = trans_from_master(0, 0, g_a);
	gauss_points_global[5] = trans_from_master(0, 0, -g_a);

	for (int i = 0; i < 6; i++)
		gauss_weights[i] = w_1;

	gauss_points_global[6] = trans_from_master(g_a, g_a, 0);
	gauss_points_global[7] = trans_from_master(g_a, -g_a, 0);
	gauss_points_global[8] = trans_from_master(-g_a, g_a, 0);
	gauss_points_global[9] = trans_from_master(-g_a, -g_a, 0);

	gauss_points_global[10] = trans_from_master(g_a, 0, g_a);
	gauss_points_global[11] = trans_from_master(g_a, 0, -g_a);
	gauss_points_global[12] = trans_from_master(-g_a, 0, g_a);
	gauss_points_global[13] = trans_from_master(-g_a, 0, -g_a);

	gauss_points_global[14] = trans_from_master(0, g_a, g_a);
	gauss_points_global[15] = trans_from_master(0, g_a, -g_a);
	gauss_points_global[16] = trans_from_master(0, -g_a, g_a);
	gauss_points_global[17] = trans_from_master(0, -g_a, -g_a);

	for (int i = 6; i < 18; i++)
		gauss_weights[i] = w_2;

	gauss_points_global[18] = trans_from_master(g_b, g_b, g_b);
	gauss_points_global[19] = trans_from_master(g_b, g_b, -g_b);
	gauss_points_global[20] = trans_from_master(g_b, -g_b, g_b);
	gauss_points_global[21] = trans_from_master(g_b, -g_b, -g_b);

	gauss_points_global[22] = trans_from_master(-g_b, g_b, g_b);
	gauss_points_global[23] = trans_from_master(-g_b, g_b, -g_b);
	gauss_points_global[24] = trans_from_master(-g_b, -g_b, g_b);
	gauss_points_global[25] = trans_from_master(-g_b, -g_b, -g_b);

	for (int i = 18; i < 26; i++)
		gauss_weights[i] = w_3;

	gauss_points_global[26] = trans_from_master(g_c, g_c, g_c);
	gauss_points_global[27] = trans_from_master(g_c, g_c, -g_c);
	gauss_points_global[28] = trans_from_master(g_c, -g_c, g_c);
	gauss_points_global[29] = trans_from_master(g_c, -g_c, -g_c);

	gauss_points_global[30] = trans_from_master(-g_c, g_c, g_c);
	gauss_points_global[31] = trans_from_master(-g_c, g_c, -g_c);
	gauss_points_global[32] = trans_from_master(-g_c, -g_c, g_c);
	gauss_points_global[33] = trans_from_master(-g_c, -g_c, -g_c);

	for (int i = 26; i < 34; i++)
		gauss_weights[i] = w_4;

	// Проверка весов
	double tmp = 0;
	for (auto w_i : gauss_weights) {
		tmp += w_i;
	}

}

point brick::get_center() {
	return point (nodes[0].x + hx / 2.0, nodes[0].y + hy / 2.0, nodes[0].z + hz / 2.0);

}

point brick::get_node(size_t i) {
	return point(nodes[i]);
}

void brick::set_env(double w_s, double sigma_s) {
	w = w_s;
	sigma = sigma_s;
}

bool brick::in_element(double x, double y, double z) {
	return
		nodes[0].x >= x && nodes[7].x <= x &&
		nodes[0].y >= y && nodes[7].y <= y &&
		nodes[0].z >= z && nodes[7].z <= z;
	
}

cmatrix(3) brick::get_matrix_value(tfunc3d G, point pn) {
	cmatrix(3) res;
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			res[i][j] = 0;
	bool inside = in_element(pn.x, pn.y, pn.z);
	array<vec3d, 3> res_in_vec;
	for(int i = 0; i < 3; i++) {

		if (inside)
			res_in_vec[i] = basis(i, pn.x, pn.y, pn.z);
		else
			res_in_vec[i] = vec3d(0, 0, 0);

		vfunc3d integ_func = [&](double x, double y, double z)->vec3d {
			cmatrix(3) G_in_point = G(x, y, z);
			vec3d func_res = G_in_point * basis(i, x, y, z);
			return func_res;
		};
		vec3d integ_part = dcomplex(0, 1) * w * sigma * integrate(integ_func);
		res_in_vec[i] = res_in_vec[i] - integ_part;

		for(int j = 0; j < 3; j++)
			res[j][i] = res_in_vec[i][j];
	}

	return res;
}

vec3d brick::basis(dof_type loc_dof, double x, double y, double z) {
	vec3d res;

	switch(loc_dof) {
		case 0: res = vec3d(1, 0, 0); break;
		case 1: res = vec3d(0, 1, 0); break;
		case 2: res = vec3d(0, 0, 1); break;
	}

	return res;
}