#include "elements_classes.h"


// ========  Отрезки ========
sector::sector() {

}
sector::sector(vector<node> nodes_s, vector<dof_type> s_dofs) {

}

double sector::L2_diff(func3d f, vector<double>& q_loc){
	return 0;
}

// ======== Прямоугольники ========

quadelement::quadelement() {
}

quadelement::quadelement(vector<node> nodes_s, vector<dof_type> s_dofs) {
	for(int i = 0; i < element_nodes; i++)
		node_array[i] = nodes_s[i];

	dofs = s_dofs;
	dofs_number = dofs.size();

	init_coords();

}

void quadelement::init_coords() {

	// Определим направление плоскости
	vec3d g1(node_array[0], node_array[1]); 
	vec3d g2(node_array[0], node_array[2]); 
	vec3d e2 = g2 - ((g1*g2) / (g1*g1)) * g1;
	vec3d e3 = g1.cross(e2);
	normal_vector = e3;

	tau[0] = g1;
	tau[1] = g2;
	tau[2] = vec3d(node_array[1], node_array[2]);

	double h1 = g1.norm(), h2 = e2.norm();

	vec3d e1 = g1 / g1.norm();
	e2 = e2 / e2.norm();
	e3 = e3 / e3.norm();

	for(int i = 0; i < 3; i++) {
		transition_matrix[0][i] = e1[i];
		transition_matrix[1][i] = e2[i];
		transition_matrix[2][i] = e3[i];
	}

	// Переведём все точки в локальные координаты
	u_0 = 0;
	v_0 = 0;
	for(int i = 0; i < 4; i++) {
		node_array_local[i] = to_local_cord(node_array[i]);
		// Найдём димаетрально противоположную точку
		if (abs(node_array_local[i].x) > GEOCONST && abs(node_array_local[i].y) > GEOCONST) {
			h_u = node_array_local[i].x;
			h_v = node_array_local[i].y;
		}
	}

	// Преобразуем точки интегрирования
	integ_point* points;
	double* weights;
	get_quadrature_rule(ELEM_TYPE::QUAD, 7, points, weights, gauss_points_n);

	gauss_points.resize(gauss_points_n);
	gauss_points_global.resize(gauss_points_n);
	gauss_weights.resize(gauss_points_n);

	// Переведём в локальные координаты
	for(int i = 0; i < gauss_points_n; i++) {
		gauss_weights[i] = weights[i];
		gauss_points[i].x = h_u * (points[i].xi + 1) / 2.0;
		gauss_points[i].y = h_v * (points[i].eta + 1) / 2.0;
		gauss_points[i].z = 0;

		gauss_points_global[i] = to_global_cord(gauss_points[i]);
	}
	
}

void quadelement::renumerate_dofs(map<dof_type, dof_type>& glob_to_loc) {
	for(int dof_i = 0; dof_i < dofs_number; dof_i++)
		dofs[dof_i] = glob_to_loc[dofs[dof_i]];
}

point quadelement::to_local_cord(point p_glob) {
	point p_shift = p_glob;

	//сдвиг

	for(int i = 0; i < 3; i++)
		p_shift[i] -= node_array[0][i];

	point p_loc;

	//поворот
	for(int i = 0; i < 3; i++) {
		p_loc[i] = 0;
		for(int j = 0; j < 3; j++)
			p_loc[i] += transition_matrix[i][j] * p_shift[j];
	}

	return p_loc;
}

point quadelement::to_global_cord(point p_loc) {
	point p_turn;

	//поворот
	for(int i = 0; i < 3; i++) {
		p_turn[i] = 0;
		for(int j = 0; j < 3; j++)
			p_turn[i] += transition_matrix[j][i] * p_loc[j];
	}

	point p_glob; 

	//сдвиг
	for(int i = 0; i < 3; i++)
		p_glob[i] = p_turn[i] + node_array[0][i];

	return p_glob;
}

dof_type& quadelement::operator [] (int i) {
	return dofs[i];
}

vector<dof_type> quadelement::get_dofs() {
	return dofs;
}

vec3d quadelement::vector_basis_v(int i, double x, double y, double z) {

	double u = get_u(x, y, z);
	double v = get_v(x, y, z);
	vec3d val;

	switch(i) {
		case 0:
			val = vec3d(0, (1-u)/2, 0);
			break;
		case 1:
			val = vec3d(0, (1+u)/2, 0);
			break;
		case 2:
			val = vec3d((1-v)/2, 0, 0);
			break;
		case 3:
			val = vec3d((1+v)/2, 0, 0);
			break;
		case 4:
			val = vec3d(0, (1-u)*v/2, 0);
			break;
		case 5:
			val = vec3d(0, (1+u)*v/2, 0);
			break;
		case 6:
			val = vec3d((1-v)*u/2, 0, 0);
			break;
		case 7:
			val = vec3d((1+v)*u/2, 0 ,0);
			break;
		case 8:
			val = vec3d(0, 1-u*u, 0);
			break;
		case 9:
			val = vec3d(0, (1-u*u)*v, 0);
			break;
		case 10:
			val = vec3d(1-v*v, 0, 0);
			break;
		case 11:
			val = vec3d((1-v*v)*u, 0, 0);
			break;

	};
	return val;

}

double quadelement::get_u(double x, double y, double z) {
	point local_p = to_local_cord(point(x, y, z));
	double u;
	u = 2 * (local_p.x - u_0) / h_u - 1;	

	return u;
}

double quadelement::get_v(double x, double y, double z) {
	point local_p = to_local_cord(point(x, y, z));
	double v;
	v = 2 * (local_p.y - v_0) / h_v - 1;	

	return v;
}


double quadelement::integrate(func3d integ_func) {
	throw;
	double res = 0;

	for(int i = 0; i < gauss_points_n; i++) {
		res += gauss_weights[i] * integ_func(gauss_points_global[i].x, gauss_points_global[i].y, gauss_points_global[i].z);
	}

	res *= jacobian;

	return res;

}

dyn_matrix quadelement::get_local_matrix(double mu) {
	throw;
	dyn_matrix M;

	M.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++) {
		M[i].resize(dofs_number);
		for(int j = 0; j <= i; j++) {
			M[i][j] = integrate([&](double x, double y, double z)->double {
				return mu * vector_basis_v(i, x, y, z) * vector_basis_v(j, x, y, z);
			});
			M[j][i] = M[i][j];
		}
	}

	return M;
}


vector<double> quadelement::get_local_right_part(vfunc3d rp_func) {
	vector<double> b;
	b.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++)
		b[i] = integrate([&](double x, double y, double z)->double {
			return  rp_func(x,y,z) * vector_basis_v(i, x, y, z);
	});

	return b;
}

cubeelement::cubeelement() {
}

cubeelement::cubeelement(vector<node> nodes_s, vector<dof_type> s_dofs) {
	for(int i = 0; i < element_nodes; i++)
		node_array[i] = nodes_s[i];

	dofs = s_dofs;
	dofs_number = dofs.size();

	init_coords();

}

double cubeelement::integrate(func3d integ_func) {
	double res = 0;

	for(int i = 0; i < gauss_points_n; i++) {
		res += gauss_weights[i] * integ_func(gauss_points_global[i].x, gauss_points_global[i].y, gauss_points_global[i].z);
	}

	res *= jacobian;

	return res;

}

void cubeelement::init_coords() {

	// Упорядочим все точки, сначала по z, потом по y, потом по x

	// по z
	sort(node_array.begin(), node_array.end(), [](node& a, node& b)->bool {	return a.z > b.z;});
	// по y
	auto y_border1 = node_array.begin()+3;
	auto y_border2 = y_border1 + 1;
	sort(node_array.begin(), y_border1, [](node& a, node& b)->bool { return a.y > b.y;});
	sort(y_border2, node_array.end(), [](node& a, node& b)->bool { return a.y > b.y;});

	// по x
	for (auto x_b = node_array.end(); x_b != node_array.end(); x_b += 2) {
		sort(x_b, x_b+1, [](node& a, node& b)->bool { return a.x > b.x;});
	}

	//Вычисляем шаги и начальную точку
	x_0 = node_array[0].x;
	y_0 = node_array[0].y;
	z_0 = node_array[0].z;

	h_x = node_array[7].x - x_0;
	h_y = node_array[7].y - y_0;
	h_z = node_array[7].z - z_0;

	// Преобразуем точки интегрирования
	integ_point* points;
	double* weights;
	get_quadrature_rule(ELEM_TYPE::CUBE, 7, points, weights, gauss_points_n);

	gauss_points.resize(gauss_points_n);
	gauss_points_global.resize(gauss_points_n);
	gauss_weights.resize(gauss_points_n);

	// Переведём в локальные координаты
	for(int i = 0; i < gauss_points_n; i++) {
		gauss_weights[i] = weights[i];
		gauss_points[i].x = points[i].xi;
		gauss_points[i].y = points[i].eta;
		gauss_points[i].z = points[i].zeta;

		gauss_points_global[i] = to_global_cord(gauss_points[i]);
	}

	jacobian = h_x * h_y * h_z / 8.0;
	
}

point cubeelement::to_global_cord(point p_loc) {
	return point(h_x * (p_loc.x + 1)/2 + x_0, h_y * (p_loc.y + 1)/2 + y_0, h_z * (p_loc.z + 1)/2 + z_0);

}

point cubeelement::to_local_cord(point p_glob) {
	return point(2 * (p_glob.x - x_0) / h_x - 1, 2 * (p_glob.y - y_0) / h_y - 1, 2 * (p_glob.z - z_0) / h_z - 1);
}

dyn_matrix cubeelement::get_local_matrix(double mu, double k_sq) {
	dyn_matrix M;

	M.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++) {
		M[i].resize(dofs_number);
		for(int j = 0; j <= i; j++) {
			M[i][j] = integrate([&](double x, double y, double z)->double {
				vec3d m1 = vector_basis_v(i, x, y, z);
				vec3d m2 = vector_basis_v(j, x, y, z);
				vec3d g1 = vector_basis_rot_v(i, x, y, z);
				vec3d g2 = vector_basis_rot_v(i, x, y, z);

				return mu * g1 * g2 - k_sq * m1 * m2;
			});
			M[j][i] = M[i][j];
		}
	}

	return M;
}


vec3d cubeelement::vector_basis_v(int i, double x, double y, double z) {
	point p_loc = to_local_cord(point(x, y, z));

	double ksi = p_loc.x;
	double etta = p_loc.y;
	double dzeta = p_loc.z;

	unordered_map<string, double> values;
	values["ksi_+"] = (1 + ksi)/2;
	values["ksi_-"] = (1 - ksi)/2;
	values["ksi_2"] = 1 - ksi*ksi;

	values["etta_+"] = (1 + etta)/2;
	values["etta_-"] = (1 - etta)/2;
	values["etta_2"] = 1 - etta*etta;

	values["dzeta_+"] = (1 + dzeta)/2;
	values["dzeta_-"] = (1 - dzeta)/2;
	values["dzeta_2"] = 1 - dzeta*dzeta;

#define phi(x, b) values[#x"_"#b] 


	vec3d val;

	switch(i) {
		// Первый порядок, первый тип
		case 0:
			val = vec3d(phi(etta, -)*phi(dzeta, -), 0, 0);
			break;
		case 1:
			val = vec3d(phi(etta, +)*phi(dzeta, -), 0, 0);
			break;
		case 2:
			val = vec3d(phi(etta, -)*phi(dzeta, +), 0, 0);
			break;
		case 3:
			val = vec3d(phi(etta, +)*phi(dzeta, +), 0, 0);
			break;
		case 4:
			val = vec3d(0, phi(ksi, -) * phi(dzeta, -), 0);
			break;
		case 5:
			val = vec3d(0, phi(ksi, +) * phi(dzeta, -), 0);
			break;
		case 6:
			val = vec3d(0, phi(ksi, -) * phi(dzeta, +), 0);
			break;
		case 7:
			val = vec3d(0, phi(ksi, +) * phi(dzeta, +), 0);
			break;
		case 8:
			val = vec3d(0, 0, phi(ksi, -) * phi(etta, -));
			break;
		case 9:
			val = vec3d(0, 0, phi(ksi, +) * phi(etta, -));
			break;
		case 10:
			val = vec3d(0, 0, phi(ksi, -) * phi(etta, +));
			break;
		case 11:
			val = vec3d(0, 0, phi(ksi, +) * phi(etta, -));
			break;

		// Первый порядок, второй тип
		case 12:
			val = ksi * vec3d(phi(etta, -)*phi(dzeta, -), 0, 0);
			break;
		case 13:
			val = ksi * vec3d(phi(etta, +)*phi(dzeta, -), 0, 0);
			break;
		case 14:
			val = ksi * vec3d(phi(etta, -)*phi(dzeta, +), 0, 0);
			break;
		case 15:
			val = ksi * vec3d(phi(etta, +)*phi(dzeta, +), 0, 0);
			break;
		case 16:
			val = etta * vec3d(0, phi(ksi, -) * phi(dzeta, -), 0);
			break;
		case 17:
			val = etta * vec3d(0, phi(ksi, +) * phi(dzeta, -), 0);
			break;
		case 18:
			val = etta * vec3d(0, phi(ksi, -) * phi(dzeta, +), 0);
			break;
		case 19:
			val = etta * vec3d(0, phi(ksi, +) * phi(dzeta, +), 0);
			break;
		case 20:
			val = dzeta * vec3d(0, 0, phi(ksi, -) * phi(etta, -));
			break;
		case 21:
			val = dzeta *vec3d(0, 0, phi(ksi, +) * phi(etta, -));
			break;
		case 22:
			val = dzeta *vec3d(0, 0, phi(ksi, -) * phi(etta, +));
			break;
		case 23:
			val = dzeta *vec3d(0, 0, phi(ksi, +) * phi(etta, -));
			break;

		// Второй порядок, первый тип
			// Грани
		case 24:
			val = vec3d(0, phi(dzeta, 2) * phi(ksi, -), 0);
			break;
		case 25:
			val = vec3d(0, 0, phi(etta, 2) * phi(ksi, -));
			break;
		case 26:
			val = vec3d(0, phi(dzeta, 2) * phi(ksi, +), 0);
			break;
		case 27:
			val = vec3d(0, 0, phi(etta, 2) * phi(ksi, +));
			break;
		case 28:
			val = vec3d(phi(dzeta, 2) * phi(etta, -), 0, 0);
			break;
		case 29:
			val = vec3d(0, 0, phi(ksi, 2) * phi(etta, -));
			break;
		case 30:
			val = vec3d(phi(dzeta, 2) * phi(etta, +), 0, 0);
			break;
		case 31:
			val = vec3d(0, 0, phi(ksi, 2) * phi(etta, +));
			break;
		case 32:
			val = vec3d(phi(etta, 2) * phi(dzeta, -), 0, 0);
			break;
		case 33:
			val = vec3d(0, phi(ksi, 2) * phi(dzeta, -), 0);
			break;
		case 34:
			val = vec3d(phi(etta, 2) * phi(dzeta, +), 0, 0);
			break;
		case 35:
			val = vec3d(0, phi(ksi, 2) * phi(dzeta, +), 0);
			break;

			//Объём
		case 36:
			val = vec3d(phi(etta, 2) * phi(dzeta, 2), 0, 0);
			break;
		case 37:
			val = vec3d(0, phi(ksi, 2) * phi(dzeta, 2), 0);
			break;
		case 38:
			val = vec3d(0, 0, phi(ksi, 2) * phi(etta, 2));
			break;

	};
#undef phi
	return val;

}

vec3d cubeelement::vector_basis_rot_v(int i, double x, double y, double z) {
	point p_loc = to_local_cord(point(x, y, z));

	double ksi = p_loc.x;
	double etta = p_loc.y;
	double dzeta = p_loc.z;

	unordered_map<string, double> values;
	values["ksi_+"] = (1 + ksi)/2;
	values["ksi_-"] = (1 - ksi)/2;
	values["ksi_2"] = 1 - ksi*ksi;

	values["ksi_+d"] =  1.0 / h_x;
	values["ksi_-d"] =  -1.0 / h_x;
	values["ksi_2d"] = -2*ksi / h_x;

	values["etta_+"] = (1 + etta)/2;
	values["etta_-"] = (1 - etta)/2;
	values["etta_2"] = 1 - etta*etta;

	values["etta_+d"] = 1.0 / h_y;
	values["etta_-d"] = -1.0 / h_y;
	values["etta_2d"] = -2 * etta / h_y;

	values["dzeta_+"] = (1 + dzeta)/2;
	values["dzeta_-"] = (1 - dzeta)/2;
	values["dzeta_2"] = 1 - dzeta*dzeta;

	values["dzeta_+d"] = 1.0 / h_z;
	values["dzeta_-d"] = -1.0 / h_z;
	values["dzeta_2d"] = -2.0 * dzeta / h_z;

#define phi(x, b) values[#x"_"#b] 
#define dphi(x, b) values[#x"_"#b"d"] 

	vec3d val;

	switch(i) {
		// Первый порядок, первый тип
		case 0:
			val = vec3d(0, phi(etta, -)*dphi(dzeta, -), -dphi(etta, -)*phi(dzeta, -));
			break;
		case 1:
			val = vec3d(0, phi(etta, +)*dphi(dzeta, -), -dphi(etta, +)*phi(dzeta, -));
			break;
		case 2:
			val = vec3d(0, phi(etta, -)*dphi(dzeta, +), -dphi(etta, -)*phi(dzeta, +));
			break;
		case 3:
			val = vec3d(0, phi(etta, +)*dphi(dzeta, +), -dphi(etta, +)*phi(dzeta, +));
			break;
		case 4:
			val = vec3d(-phi(ksi, -) * dphi(dzeta, -), 0, dphi(ksi, -) * dphi(dzeta, -));
			break;
		case 5:
			val = vec3d(-phi(ksi, +) * dphi(dzeta, -), 0, dphi(ksi, +) * dphi(dzeta, -));
			break;
		case 6:
			val = vec3d(-phi(ksi, -) * dphi(dzeta, +), 0, dphi(ksi, -) * dphi(dzeta, +));
			break;
		case 7:
			val = vec3d(-phi(ksi, +) * dphi(dzeta, +), 0, dphi(ksi, +) * dphi(dzeta, +));
			break;
		case 8:
			val = vec3d(phi(ksi, -) * dphi(etta, -), -dphi(ksi, -) * phi(etta, -), 0);
			break;
		case 9:
			val = vec3d(phi(ksi, +) * dphi(etta, -), -dphi(ksi, +) * phi(etta, -), 0);
			break;
		case 10:
			val = vec3d(phi(ksi, -) * dphi(etta, +), -dphi(ksi, -) * phi(etta, +), 0);
			break;
		case 11:
			val = vec3d(phi(ksi, +) * dphi(etta, +), -dphi(ksi, +) * phi(etta, +), 0);
			break;

		// Первый порядок, второй тип
		case 12:
			val = ksi * vec3d(0, phi(etta, -)*dphi(dzeta, -), -dphi(etta, -)*phi(dzeta, -));
			break;
		case 13:
			val = ksi * vec3d(0, phi(etta, +)*dphi(dzeta, -), -dphi(etta, +)*phi(dzeta, -));
			break;
		case 14:
			val = ksi * vec3d(0, phi(etta, -)*dphi(dzeta, +), -dphi(etta, -)*phi(dzeta, +));
			break;
		case 15:
			val = ksi * vec3d(0, phi(etta, +)*dphi(dzeta, +), -dphi(etta, +)*phi(dzeta, +));
			break;
		case 16:
			val = etta * vec3d(-phi(ksi, -) * dphi(dzeta, -), 0, dphi(ksi, -) * dphi(dzeta, -));
			break;
		case 17:
			val = etta * vec3d(-phi(ksi, +) * dphi(dzeta, -), 0, dphi(ksi, +) * dphi(dzeta, -));
			break;
		case 18:
			val = etta * vec3d(-phi(ksi, -) * dphi(dzeta, +), 0, dphi(ksi, -) * dphi(dzeta, +));
			break;
		case 19:
			val = etta * vec3d(-phi(ksi, +) * dphi(dzeta, +), 0, dphi(ksi, +) * dphi(dzeta, +));
			break;
		case 20:
			val = dzeta * vec3d(phi(ksi, -) * dphi(etta, -), -dphi(ksi, -) * phi(etta, -), 0);
			break;
		case 21:
			val = dzeta * vec3d(phi(ksi, +) * dphi(etta, -), -dphi(ksi, +) * phi(etta, -), 0);
			break;
		case 22:
			val = dzeta * vec3d(phi(ksi, -) * dphi(etta, +), -dphi(ksi, -) * phi(etta, +), 0);
			break;
		case 23:
			val = dzeta * vec3d(phi(ksi, +) * dphi(etta, +), -dphi(ksi, +) * phi(etta, +), 0);
			break;

		// Второй порядок, первый тип
			// Грани
		case 24:
			val = vec3d(-dphi(dzeta, 2) * phi(ksi, -), 0, phi(dzeta, 2) * dphi(ksi, -));
			break;
		case 25:
			val = vec3d(dphi(etta, 2) * phi(ksi, -), - phi(etta, 2,) * dphi(ksi, -), 0);
			break;
		case 26:
			val = vec3d(-dphi(dzeta, 2) * phi(ksi, +), 0, phi(dzeta, 2) * dphi(ksi, +));
			break;
		case 27:
			val = vec3d(dphi(etta, 2) * phi(ksi, +), - phi(etta, 2,) * dphi(ksi, +), 0);
			break;
		case 28:
			val = vec3d(0, dphi(dzeta, 2) * phi(etta, -), -phi(dzeta, 2) * dphi(etta, -));
			break;
		case 29:
			val = vec3d(phi(ksi, 2) * dphi(etta, -), -dphi(ksi, 2) * phi(etta, -), 0);
			break;
		case 30:
			val = vec3d(0, dphi(dzeta, 2) * phi(etta, +), -phi(dzeta, 2) * dphi(etta, +));
			break;
		case 31:
			val = vec3d(phi(ksi, 2) * dphi(etta, +), -dphi(ksi, 2) * phi(etta, +), 0);
			break;
		case 32:
			val = vec3d(0, phi(etta, 2) * dphi(dzeta, -), -dphi(etta, 2) * phi(dzeta, -));
			break;
		case 33:
			val = vec3d(-phi(ksi, 2) * dphi(dzeta, -), 0, dphi(ksi, 2) * phi(dzeta, -));
			break;
		case 34:
			val = vec3d(0, phi(etta, 2) * dphi(dzeta, +), -dphi(etta, 2) * phi(dzeta, +));
			break;
		case 35:
			val = vec3d(-phi(ksi, 2) * dphi(dzeta, +), 0, dphi(ksi, 2) * phi(dzeta, +));
			break;

			//Объём
		case 36:
			val = vec3d(0, phi(etta, 2) * dphi(dzeta, 2), -dphi(etta, 2) * phi(dzeta, 2));
			break;
		case 37:
			val = vec3d(-phi(ksi, 2) * dphi(dzeta, 2), 0, dphi(ksi, 2) * phi(dzeta, 2));
			break;
		case 38:
			val = vec3d(phi(ksi, 2) * dphi(etta, 2), -dphi(ksi, 2) * phi(etta, 2), 0);
			break;

	};
#undef phi
	return val;

}

/*double quadelement::L2_diff(func3d f, vector<double>& q_loc){

	int loc_n = q_loc.size();

	double res = integrate([&](double x, double y, double z){
		double u = 0;
		for (int i = 0; i < loc_n; i++) {
			u += q_loc[i] * (this->*scalar_basis[i])(x, y, z);
		}

		double f_v = f(x, y, z);
		return (u - f_v)*(u - f_v);
	});

	return 0;
}

double trelement::vector_jump_L2(vfunc3d f1, vfunc3d f2) {
	return integrate([&](double x, double y, double z)->double {
		vec3d f1_v = f1(x,y,z);
		vec3d f2_v = f2(x,y,z);
		double diff = (f1_v - f2_v)*normal_vector;

		return diff*diff;
	});
}
*/


// ======== Тетраэдры ========
/*
tetelement::tetelement() {
}

tetelement::tetelement(vector<node> nodes_s, vector<dof_type> s_dofs) {
	node_array[0] = nodes_s[0];
	node_array[1] = nodes_s[1];
	node_array[2] = nodes_s[2];
	node_array[3] = nodes_s[3];

	dofs = s_dofs;
	dofs_number = dofs.size();

	init_coords();
}

int& tetelement::operator [] (int i) {
	return edge_array[i];
}


int tetelement::get_ph_area() {
	return ph_area;
}

void tetelement::set_ph_area(int sph_area) {
	ph_area = sph_area;
}

node tetelement::get_local_node(int i) {
	return node_array[i];
}

point tetelement::get_center() {
	point cent(0, 0, 0);
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 4; j++)
			cent[i] += node_array[j][i] / 4.0;

	return cent;

}

double tetelement::scalar_basis_v(int i, double x, double y, double z) {
	return (this->*scalar_basis[i])(x,y,z); 
}

vec3d tetelement::scalar_basis_grad_v(int i, double x, double y, double z) {
	return (this->*scalar_basis_grad[i])(x,y,z); 
}

vector<dof_type> tetelement::get_dofs() {
	return dofs;
}


bool tetelement::in_element(double x, double y, double z) {
	point p_glob(x,y,z);

#ifdef DEBUGOUTP
	array<double, 4> lambdas;
	for(int i = 0; i < 4; i++) {
		lambdas[i] = lambda(i, p_glob);
	}
	for(int i = 0; i < 4; i++) {
		if(lambdas[i] > 1 + GEOCONST|| lambdas[i] < 0 - GEOCONST)
			return false;
	}
#else
	for(int i = 0; i < 4; i++) {
		double L = lambda(i, p_glob);
		if(L > 1 + GEOCONST|| L < 0 - GEOCONST)
			return false;
	}
#endif

	return true;

}
bool tetelement::valid_for_tree_node(double x0, double x1, double y0, double y1, double z0, double z1) {

	//Проерка на характеристические точки в том числе, если тетраэдр лежит внутри куба (ну вытянутого)
	for(int i = 0; i < 5; i++)
		if(ch_points[i][0] >= x0 && ch_points[i][0] <= x1 && ch_points[i][1] >= y0 && ch_points[i][1] <= y1 && ch_points[i][2] >= z0 && ch_points[i][2] <= z1)
			return true;

	//Куб лежит внутри тэтраэдра (хотя бы частично)
	if(in_element(x0, y0, z0) || in_element(x1, y0, z0) || in_element(x0, y1, z0) || in_element(x1, y1, z0) ||
		in_element(x0, y0, z1) || in_element(x1, y0, z1) || in_element(x0, y1, z1) || in_element(x1, y1, z1))
		return true;

	//Самое весёлое - куб пересекает тэтраэдр
	double t[6][6]; //параметры прямых для всех прямых (6) и всех плоскостей (6)
	double cords_t[6][6][3]; //прямая, плоскость, координата

	for(int i = 0; i < 6; i++) {
		t[i][0] = (x0 - edges_b[i][0]) / edges_a[i][0]; // x = x0
		t[i][1] = (x1 - edges_b[i][0]) / edges_a[i][0]; // x = x1

		t[i][2] = (y0 - edges_b[i][1]) / edges_a[i][1]; // y = y0
		t[i][3] = (y1 - edges_b[i][1]) / edges_a[i][1]; // y = y1

		t[i][4] = (z0 - edges_b[i][2]) / edges_a[i][2]; // z = z0
		t[i][5] = (z1 - edges_b[i][2]) / edges_a[i][2]; // z = z1

	}
	
	for(int i = 0; i < 6; i++)
		for(int j = 0; j < 6; j++) 
			for(int k = 0; k < 3; k++) 
				cords_t[i][j][k] = edges_a[i][k] * t[i][j] + edges_b[i][k];

	for(int i = 0; i < 6; i++) { //берём прямоую и проверяем, что пересечение с плоскостью попадает в раcсматриваемый отрезок прямой и в рассатриваемую часть плоскости

		if(	   t[i][0] >= 0 && t[i][0] <= 1 && cords_t[i][0][1] >= y0 && cords_t[i][0][1] >= y1 && cords_t[i][0][2] >= z0 && cords_t[i][0][2] <= z1	// x = x0
			|| t[i][1] >= 0 && t[i][1] <= 1 && cords_t[i][1][1] >= y0 && cords_t[i][1][1] >= y1 && cords_t[i][1][2] >= z0 && cords_t[i][1][2] <= z1 // x = x1
			|| t[i][2] >= 0 && t[i][2] <= 1 && cords_t[i][2][0] >= x0 && cords_t[i][2][0] <= x1 && cords_t[i][2][2] >= z0 && cords_t[i][2][2] <= z1 // y = y0
			|| t[i][3] >= 0 && t[i][3] <= 1 && cords_t[i][3][0] >= x0 && cords_t[i][3][0] <= x1 && cords_t[i][3][2] >= z0 && cords_t[i][3][2] <= z1 // y = y1
			|| t[i][4] >= 0 && t[i][4] <= 1 && cords_t[i][4][0] >= x0 && cords_t[i][4][0] <= x1 && cords_t[i][4][1] >= y0 && cords_t[i][4][1] <= y1 // z = z0
			|| t[i][5] >= 0 && t[i][5] <= 1 && cords_t[i][5][0] >= x0 && cords_t[i][5][0] <= x1 && cords_t[i][5][1] >= y0 && cords_t[i][5][1] <= y1 // z = z1
			)
			return true;


	}


	return false;

}

void tetelement::init_coords() {

	

	//расчёт L-координат
	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 3; j++)
			D_matrix[j][i] = node_array[i][j];

	for(int i = 0; i < 4; i++)
		D_matrix[3][i] = 1;

	L_cord_matrix = inverse4(D_matrix, det_D);

	jacobian = fabs(det_D);

	//Точки Гаусса на мастер-элементе
	double gauss_a = (5.0 - sqrt(5.0)) / 20.0;
	double gauss_b = (5.0 + 3.0*sqrt(5.0)) / 20.0;

	double Gauss_cord[4][4];
	double Gauss_cord_gl[4][4];

	//i-й столбец, i-я точка

	Gauss_cord[0][0] = 1 - gauss_b - 2*gauss_a;
	Gauss_cord[1][0] = gauss_b;
	Gauss_cord[2][0] = Gauss_cord[3][0] = gauss_a;

	Gauss_cord[0][1] = 1 - gauss_b - 2*gauss_a;
	Gauss_cord[1][1] = Gauss_cord[3][1] = gauss_a;
	Gauss_cord[2][1] = gauss_b;

	Gauss_cord[0][2] = 1 - gauss_b - 2*gauss_a;
	Gauss_cord[1][2] = Gauss_cord[2][2] = gauss_a;
	Gauss_cord[3][2] = gauss_b;

	Gauss_cord[0][3] = 1 - 3*gauss_a;
	Gauss_cord[1][3] = Gauss_cord[2][3] = Gauss_cord[3][3] = gauss_a;

	//Перевод на текущий тетраэдр

	for(int i = 0; i < 4; i++) {
		for(int j = 0 ; j < gauss_points_tet; j++) {
			Gauss_cord_gl[i][j] = 0;
			for(int k = 0; k < 4; k++)
				Gauss_cord_gl[i][j] += D_matrix[i][k] * Gauss_cord[k][j];
		}
	}

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 3; j++)
			gauss_points[i][j] = Gauss_cord_gl[j][i];

	for(int i = 0; i < gauss_points_tet; i++)
		gauss_weights[i] = 1.0 / 24.0;

	//Ниже - для данные для построения дерева поиска

	for(int i = 0; i < 3; i++) {
		ch_points[0][i] = 0;
		for(int j = 0; j < 4; j++)
			ch_points[0][i] += node_array[j][i] / 4.0;
	}

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 3; j++)
			ch_points[i+1][j] = node_array[i][j];

	//Представление прямых тетраэдра в параметричком виде a*t + b, сам отрезок - ребо получается при 0<=t<=1
	//Нужно для дерева поиска

	for(int i = 0; i < 3; i++) {
		edges_a[0][i] = node_array[1][i] - node_array[0][i];
		edges_b[0][i] = node_array[0][i];

		edges_a[1][i] = node_array[2][i] - node_array[0][i];
		edges_b[1][i] = node_array[0][i];

		edges_a[2][i] = node_array[3][i] - node_array[0][i];
		edges_b[2][i] = node_array[0][i];

		edges_a[3][i] = node_array[2][i] - node_array[1][i];
		edges_b[3][i] = node_array[1][i];

		edges_a[4][i] = node_array[3][i] - node_array[1][i];
		edges_b[4][i] = node_array[1][i];

		edges_a[5][i] = node_array[3][i] - node_array[2][i];
		edges_b[5][i] = node_array[2][i];

	}

	
	//Формирование базиса в виде массива функций
	scalar_basis.push_back(&tetelement::basis_1);
	scalar_basis.push_back(&tetelement::basis_2);
	scalar_basis.push_back(&tetelement::basis_3);
	scalar_basis.push_back(&tetelement::basis_4);

	scalar_basis_grad.push_back(&tetelement::grad_basis_1);
	scalar_basis_grad.push_back(&tetelement::grad_basis_2);
	scalar_basis_grad.push_back(&tetelement::grad_basis_3);
	scalar_basis_grad.push_back(&tetelement::grad_basis_4);

	scalar_basis.resize(dofs_number);
	scalar_basis_grad.resize(dofs_number);

}

double tetelement::lambda(int i, point p_glob) {

	double res = L_cord_matrix[i][3];
	for(int j = 0; j < 3; j++)
		res += L_cord_matrix[i][j] * p_glob[j];

	return res;

}

vec3d tetelement::grad_lambda(int i) {
	return vec3d(L_cord_matrix[i][0], L_cord_matrix[i][1], L_cord_matrix[i][2]);
}

double tetelement::integrate(func3d integ_func) {

	double res = 0;
	for(int i = 0; i < gauss_points_tet; i++) {
		double func_v = integ_func(gauss_points[i].x, gauss_points[i].y, gauss_points[i].z);
		res += gauss_weights[i] * func_v;
	}

	res *= jacobian;

	return res;

}

dyn_matrix tetelement::get_local_matrix(double mu) {
	dyn_matrix A_loc;
	A_loc.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++) {
		A_loc[i].resize(dofs_number);
		for(int j = 0; j <= i; j++) {
			A_loc[i][j] += integrate([&](double x, double y, double z)->double {
				return mu * (this->*scalar_basis_grad[i])(x,y,z) * (this->*scalar_basis_grad[j])(x,y,z);
			});
			A_loc[j][i] = A_loc[i][j];
		}
	}

	return A_loc; 
}


vector<double> tetelement::get_local_right_part(func3d rp_func) {
	vector<double> b_loc;
	b_loc.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++)
		b_loc[i] = integrate([&](double x, double y, double z)->double{
			return rp_func(x,y,z) * (this->*scalar_basis[i])(x,y,z);
	});

	return b_loc;
}

double tetelement::L2_diff(func3d f, vector<double>& q_loc, vector<double>& q_virtual){

	int loc_n = q_loc.size();
	int virt_n = q_virtual.size();

	double res = integrate([&](double x, double y, double z){
		double u = 0;
		for (int j = 0; j < loc_n; j++) {
			double basis_v = (this->*scalar_basis[j])(x, y, z);
			double q_v = 0;
			for (int i = 0; i < virt_n; i++) {
				q_v += q_virtual[i] * q_loc[j];
			}
			u += q_v * basis_v;
		}

		double f_v = f(x, y, z);
		return (u - f_v)*(u - f_v);
	});

	return res;
}


// == Базисные функции ==

double tetelement::basis_1(double x, double y, double z) {
	return lambda(0, point(x,y,z));
}

double tetelement::basis_2(double x, double y, double z) {
	return lambda(1, point(x,y,z));
}

double tetelement::basis_3(double x, double y, double z) {
	return lambda(2, point(x,y,z));
}

double tetelement::basis_4(double x, double y, double z) {
	return lambda(3, point(x,y,z));
}

vec3d tetelement::grad_basis_1(double x, double y, double z) {
	return grad_lambda(0);
}

vec3d tetelement::grad_basis_2(double x, double y, double z) {
	return grad_lambda(1);
}

vec3d tetelement::grad_basis_3(double x, double y, double z) {
	return grad_lambda(2);
}

vec3d tetelement::grad_basis_4(double x, double y, double z) {
	return grad_lambda(3);
}*/