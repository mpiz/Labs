#include "elements_classes.h"
#ifdef DEBUGOUTP
#include <iostream>
#endif


const double mu0 = 16 * atan(1.0) * 1e-7;
const double eps0 = 8.85418782 * 1e-12;


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

	jacobian = h_u * h_v / 4.0;

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

vec3d quadelement::mull_transT(vec3d v) {
	vec3d v_turn(0, 0, 0);
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++)
			v_turn[i] += transition_matrix[j][i] * v[j];
	}

	return v_turn;
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
			val = vec3d(1-v*v, 0, 0);
			break;

	};
	vec3d trans = mull_transT(val);

	return trans;

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
	double res = 0;

	for(int i = 0; i < gauss_points_n; i++) {
		res += gauss_weights[i] * integ_func(gauss_points_global[i].x, gauss_points_global[i].y, gauss_points_global[i].z);
	}

	res *= jacobian;

	return res;

}

dyn_matrix quadelement::get_local_matrix(double mu) {
	dyn_matrix M;

	M.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++) {
		M[i].resize(dofs_number);
		for(int j = 0; j <= i; j++) {
			M[i][j] = integrate([&](double x, double y, double z)->double {
				vec3d v1 = vector_basis_v(i, x, y, z);
				vec3d v2 = vector_basis_v(j, x, y, z);
				return v1 * v2 ;
			});
			M[j][i] = M[i][j];
		}
	}

	return M;
}


vector<double> quadelement::get_local_right_part(vfunc3d rp_func) {
	vector<double> b;
	b.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++) {
		b[i] = integrate([&](double x, double y, double z)->double {
			vec3d v = vector_basis_v(i, x, y, z);
			vec3d r = rp_func(x,y,z);
			double rs = r * v;
			return  r * v;
		});
	}

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

	/*// по z
	sort(node_array.begin(), node_array.end(), [](node& a, node& b)->bool {	return a.z > b.z;});
	// по y
	auto y_border1 = node_array.begin()+3;
	auto y_border2 = y_border1 + 1;
	sort(node_array.begin(), y_border1, [](node& a, node& b)->bool { return a.y > b.y;});
	sort(y_border2, node_array.end(), [](node& a, node& b)->bool { return a.y > b.y;});

	// по x
	for (auto x_b = node_array.end(); x_b != node_array.end(); x_b += 2) {
		sort(x_b, x_b+1, [](node& a, node& b)->bool { return a.x > b.x;});
	}*/

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

vector<dof_type> cubeelement::get_dofs() {
	return dofs;
}

dyn_matrix cubeelement::get_local_matrix(double mu) {
	dyn_matrix M;

	M.resize(dofs_number);
	double k_sq = -1.8e6 * eps0;
	mu /= mu0;

	for(int i = 0; i < dofs_number; i++) {
		M[i].resize(dofs_number);

		for(int j = 0; j <= i; j++) {
			M[i][j] = integrate([&](double x, double y, double z)->double {
				vec3d m1 = vector_basis_v(i, x, y, z);
				vec3d m2 = vector_basis_v(j, x, y, z);
				vec3d g1 = vector_basis_rot_v(i, x, y, z);
				vec3d g2 = vector_basis_rot_v(j, x, y, z);

				return mu * g1 * g2 - k_sq * m1 * m2;
			});
			M[j][i] = M[i][j];
		}
	}

	return M;
}

vector<double> cubeelement::get_local_right_part(vfunc3d rp_func) {
	vector<double> b;
	b.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++)
		b[i] = integrate([&](double x, double y, double z)->double {
			return  rp_func(x,y,z) * vector_basis_v(i, x, y, z);
	});

	return b;
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
	values["ksi_2d"] = -4.0*ksi / h_x;

	values["etta_+"] = (1 + etta)/2;
	values["etta_-"] = (1 - etta)/2;
	values["etta_2"] = 1 - etta*etta;

	values["etta_+d"] = 1.0 / h_y;
	values["etta_-d"] = -1.0 / h_y;
	values["etta_2d"] = -4.0 * etta / h_y;

	values["dzeta_+"] = (1 + dzeta)/2;
	values["dzeta_-"] = (1 - dzeta)/2;
	values["dzeta_2"] = 1 - dzeta*dzeta;

	values["dzeta_+d"] = 1.0 / h_z;
	values["dzeta_-d"] = -1.0 / h_z;
	values["dzeta_2d"] = -4.0 * dzeta / h_z;

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
			val = vec3d(-phi(ksi, -) * dphi(dzeta, -), 0, dphi(ksi, -) * phi(dzeta, -));
			break;
		case 5:
			val = vec3d(-phi(ksi, +) * dphi(dzeta, -), 0, dphi(ksi, +) * phi(dzeta, -));
			break;
		case 6:
			val = vec3d(-phi(ksi, -) * dphi(dzeta, +), 0, dphi(ksi, -) * phi(dzeta, +));
			break;
		case 7:
			val = vec3d(-phi(ksi, +) * dphi(dzeta, +), 0, dphi(ksi, +) * phi(dzeta, +));
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
			val = etta * vec3d(-phi(ksi, -) * dphi(dzeta, -), 0, dphi(ksi, -) * phi(dzeta, -));
			break;
		case 17:
			val = etta * vec3d(-phi(ksi, +) * dphi(dzeta, -), 0, dphi(ksi, +) * phi(dzeta, -));
			break;
		case 18:
			val = etta * vec3d(-phi(ksi, -) * dphi(dzeta, +), 0, dphi(ksi, -) * phi(dzeta, +));
			break;
		case 19:
			val = etta * vec3d(-phi(ksi, +) * dphi(dzeta, +), 0, dphi(ksi, +) * phi(dzeta, +));
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

bool cubeelement::in_element(point pn) {
	return
		pn.x >= x_0 && pn.x <= node_array[7].x &&
		pn.y >= y_0 && pn.x <= node_array[7].y &&
		pn.z >= x_0 && pn.z <= node_array[7].z;
	
}