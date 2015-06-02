#include "elements_classes.h"

const int edge_lambdas[6][2] = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}};
const int face_lambdas[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

double kernel(int order, double x) {
	double result;

	switch(order) {
	case 0 : 
		result = -2.0 * sqrt(3.0 / 2.0);
		break;
	case 1 :
		result = -2.0 * sqrt(5.0 / 2.0) * x;
		break;

	}

	return result;
}

double kernel_d(int order, double x) {
	double result;

	switch(order) {
	case 0 : 
		result = 0;
		break;
	case 1 :
		result = -2.0 * sqrt(5.0 / 2.0);
		break;

	}

	return result;
}


// ========  Отрезки ========
sector::sector() {

}
sector::sector(vector<node> nodes_s, vector<dof_type> s_dofs) {

}

double sector::L2_diff(func3d f, vector<double>& q_loc){
	return 0;
}

// ======== Треугольники ========

trelement::trelement() {
}

trelement::trelement(vector<node>& nodes_s) {
	node_array[0] = nodes_s[0];
	node_array[1] = nodes_s[1];
	node_array[2] = nodes_s[2];
}

trelement::trelement(vector<node> nodes_s, vector<dof_type> s_dofs) {
	node_array[0] = nodes_s[0];
	node_array[1] = nodes_s[1];
	node_array[2] = nodes_s[2];

	dofs = s_dofs;
	dofs_number = dofs.size();

	init_cords();
}

int& trelement::operator [] (int i) {
	return edge_array[i];
}

node trelement::local_node(int i) {
	return node_array[i];
}

int trelement::get_ph_area() {
	return ph_area;
}

void trelement::set_ph_area(int sph_area) {
	ph_area = sph_area;
}

vector<dof_type> trelement::get_dofs() {
	return dofs;
}

point trelement::get_center() {
	return center;
}

void trelement::init_cords() {

	// Определим барицентр треугольника
	for(int i = 0; i < 3; i++) {
		center[i] = 0;
		for(int j = 0; j < 3; j++)
			center[i] += node_array[j][i] / 3.0;
	}

	//Построение локальной системы координат
	vec3d g1(node_array[0], node_array[1]); 
	vec3d g2(node_array[0], node_array[2]); 
	vec3d e2 = g2 - ((g1*g2) / (g1*g1)) * g1;
	vec3d e3 = g1.cross(e2);

	tau[0] = g1;
	tau[1] = g2;
	tau[2] = vec3d(node_array[1], node_array[2]);

	double h1 = g1.norm(), h2 = e2.norm();

	trpoint[0] = point(0, 0, 0);
	trpoint[1] = point(g1.norm(), 0, 0);

	vec3d e1 = g1 / g1.norm();
	e2 = e2 / e2.norm();
	e3 = e3 / e3.norm();
	normal_vector = e3;

	for(int i = 0; i < 3; i++) {
		transition_matrix[0][i] = e1[i];
		transition_matrix[1][i] = e2[i];
		transition_matrix[2][i] = e3[i];
	}

	trpoint[2] = (transition_matrix * g2).to_point();

	generate_L_cords();



	//Нахождене точек интегрирования по Гауссу, в локальной системе координат

	jacobian = (g1.cross(g2)).norm();

	double Gauss_cord_gl[3][gauss_points_tr];

	//Перевод на текущий треугольник

	for(int i = 0; i < 3; i++) {
		for(int j = 0 ; j < gauss_points_tr; j++) {
			Gauss_cord_gl[i][j] = 0;
			for(int k = 0; k < 3; k++)
				Gauss_cord_gl[i][j] += D_matrix[i][k] * tr_integration::gauss_points_master[j][k];
		}
	}

	for(int i = 0; i < gauss_points_tr; i++)
		gauss_points[i] = point(Gauss_cord_gl[0][i], Gauss_cord_gl[1][i], 0);

	double wt  = 0;

	for(int i = 0; i < gauss_points_tr; i++) {
		gauss_points_global[i] = to_global_cord(gauss_points[i]);
		gauss_weights[i] = tr_integration::gauss_weights[i];
		wt += gauss_weights[i];
	}


}

bool trelement::in_element(double x, double y, double z) {

	point p_loc = to_local_cord(point(x,y,z));
	double lambda_v;

	// Если лежит вне плоскости, то не в элементе
	if (fabs(p_loc.z) > GEOCONST)
		return false;

	for(int i = 0; i < 3; i++) {
		lambda_v = lambda(i, p_loc);
		if(lambda_v < -GEOCONST || lambda_v > 1 + GEOCONST)
			return false;
	}

	return true;

}

point trelement::to_local_cord(point p_glob) {
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

point trelement::to_global_cord(point p_loc) {
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

vec3d trelement::to_global_cord(vec3d v_loc) {
	vec3d v_glob;

	//поворот
	for(int i = 0; i < 3; i++) {
		v_glob[i] = 0;
		for(int j = 0; j < 3; j++)
			v_glob[i] += transition_matrix[j][i] * v_loc[j];
	}


	return v_glob;

}

double trelement::scalar_basis_v(int i, double x, double y, double z) {
	point p_glob(x,y,z);
//	point p_loc = to_local_cord(p_glob);

	return (this->*scalar_basis[i])(x,y,z);

}

vec3d trelement::get_tau(int i) {
	return tau[i];
}

double trelement::integrate(func3d integ_func) {

	double res = 0;

	for(int i = 0; i < gauss_points_tr; i++) {
		res += gauss_weights[i] * integ_func(gauss_points_global[i].x, gauss_points_global[i].y, gauss_points_global[i].z);
	}

	res *= jacobian;

	return res;

}

void trelement::generate_L_cords() {

	//Формирование L-координат
	for(int i = 0; i < 2; i++)
		for(int j = 0; j < 3; j++)
			D_matrix[i][j] = trpoint[j][i];

	for(int i = 0; i < 3; i++)
		D_matrix[2][i] = 1;

	L_cord_matrix = inverse3(D_matrix, det_D);

}

double trelement::lambda(int l_i, point p_loc) {
	double res = 0;

	for(int i = 0; i < 2; i++)
		res += L_cord_matrix[l_i][i] * p_loc[i];
	res += L_cord_matrix[l_i][2];

	return res;

}

dyn_matrix trelement::get_local_matrix(double mu) {
	dyn_matrix M;

	M.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++) {
		M[i].resize(dofs_number);
		for(int j = 0; j <= i; j++) {
			M[i][j] = integrate([&](double x, double y, double z)->double {
				return mu * (this->*scalar_basis_grad[i])(x,y,z) * (this->*scalar_basis_grad[j])(x,y,z);
			});
			M[j][i] = M[i][j];
		}
	}

	return M;
}


vector<double> trelement::get_local_right_part(func3d rp_func) {
	vector<double> b;
	b.resize(dofs_number);

	for(int i = 0; i < dofs_number; i++)
		b[i] = integrate([&](double x, double y, double z)->double {
			return  rp_func(x,y,z) * (this->*scalar_basis[i])(x,y,z);
	});

	return b;
}

double trelement::L2_diff(func3d f, vector<double>& q_loc){

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


vec3d trelement::grad_lambda(int i) {
	vec3d grad;
	for(int j = 0; j < 2; j++) {
		grad[j] = L_cord_matrix[i][j];
	}

	return to_global_cord(grad);
}

// == Базисные функции ==
double trelement::basis_1(double x, double y, double z) {
	return lambda(0, to_local_cord(point(x,y,z)));
}

double trelement::basis_2(double x, double y, double z) {
	return lambda(1, to_local_cord(point(x,y,z)));
}

double trelement::basis_3(double x, double y, double z) {
	return lambda(2, to_local_cord(point(x,y,z)));
}


vec3d trelement::grad_basis_1(double x, double y, double z) {
	return grad_lambda(0);
}

vec3d trelement::grad_basis_2(double x, double y, double z) {
	return grad_lambda(1);
}

vec3d trelement::grad_basis_3(double x, double y, double z) {
	return grad_lambda(2);
}

// ======== Треугольники DG ========


trface::trface(vector<node> nodes_s) : trelement(nodes_s) {
	dofs_number = 0;
	el_count = 0;
	init_cords();
}

void trface::add_element(tetelement* el) {

	face_elements[el_count] = el;
	if(el_count == 0) {
		// Нормаль должна быть венешней для первого элемента. Убедимся в этом
		vec3d cent_vector(face_elements[0]->get_center(), center);
		double vec_dp = cent_vector * normal_vector;
		// Если надо - развернём
		if (vec_dp < 0)
			normal_vector = (-1) * normal_vector;
	}
	auto el_dofs = face_elements[el_count]->get_dofs();
	auto el_dofs_n = el_dofs.size();
	elements_dofs[el_count] = el_dofs_n;
	dofs_number += el_dofs_n;

	for(auto dof_iter = el_dofs.begin(); dof_iter != el_dofs.end(); dof_iter++)
		dofs.push_back(*dof_iter);

	el_count++;

}

vector<double> trface::get_local_right_part(func3d rp_func) {
	vector<double> b_loc;
	b_loc.resize(dofs_number);
	for(int i = 0; i < dofs_number; i++) {
		b_loc[i] = integrate([&](double x, double y, double z)->double {
			return face_elements[0]->scalar_basis_grad_v(i, x, y, z) * normal_vector * rp_func(x, y, z);
		});

	}
	return b_loc;
}

int trface::get_el_number() {
	return el_count;
}

dyn_matrix trface::get_local_matrix(double lambda) {
	int get_dofs;

	dyn_matrix A_loc;
	A_loc.resize(dofs_number);
	for(int i = 0; i < dofs_number; i++)
		A_loc[i].resize(dofs_number);

	// Диагональные блоки
	int count_dof = 0;
	for(int el_i = 0; el_i < face_el_n; el_i++) {
		for(int dof_i = 0; dof_i < elements_dofs[el_i]; dof_i++)  {
			for(int dof_j = 0; dof_j < elements_dofs[el_i]; dof_j++) {
				A_loc[count_dof+dof_i][count_dof+dof_j] = integrate([&](double x, double y, double z)->double {
					double phi1, phi2;
					vec3d v1, v2;
					phi1 = face_elements[el_i]->scalar_basis_v(dof_i, x, y, z);
					phi2 = face_elements[el_i]->scalar_basis_v(dof_j, x, y, z);
					v1 = face_elements[el_i]->scalar_basis_grad_v(dof_i, x, y, z);
					v2 = face_elements[el_i]->scalar_basis_grad_v(dof_j, x, y, z);
					vec3d nv = normal_vector;

					double res = lambda * (phi1 * v2 - phi2 * v1) * normal_vector / 2.0;
					if (el_i == 1)
						res *= -1;

					return res;

				});

			}
		}
		count_dof += elements_dofs[el_i];

	}

	// Внедиагональные блоки
	for(int dof2_i = 0; dof2_i < elements_dofs[1]; dof2_i++) {
		for(int dof1_j = 0; dof1_j < elements_dofs[0]; dof1_j++) {
			A_loc[elements_dofs[0]+dof2_i][dof1_j] = integrate([&](double x, double y, double z)->double {
					double phi1, phi2;
					vec3d v1, v2;
					phi1 = face_elements[0]->scalar_basis_v(dof1_j, x, y, z);
					phi2 = face_elements[1]->scalar_basis_v(dof2_i, x, y, z);
					v1 = face_elements[0]->scalar_basis_grad_v(dof1_j, x, y, z);
					v2 = face_elements[1]->scalar_basis_grad_v(dof2_i, x, y, z);

					double res = lambda * (phi1 * v2 + phi2 * v2) * normal_vector / 2.0;

					return res;
			});
			A_loc[dof1_j][elements_dofs[0]+dof2_i] = -A_loc[elements_dofs[0]+dof2_i][dof1_j];
		}
	}

	return A_loc;
}

// ======== Тетраэдры ========

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
	return center;
}

double tetelement::scalar_basis_v(int i, double x, double y, double z) {
	point p_glob(x, y, z);
	double result;



	// То, что заместо первого порядка
	if (i <= 3) {
		double lambda_i = lambda(i, p_glob);
		result = lambda_i * (lambda_i - 0.5);
	}
	// Рёберные функции второго порядка
	else if(i <= 9) {
		int shift = i - 3;
		double lambda_1 = lambda(edge_lambdas[shift][0], p_glob);
		double lambda_2 = lambda(edge_lambdas[shift][1], p_glob);
		result = 2 * lambda_1 * lambda_2;
	}
	// Рёберные функции третьего порядка
	else if (i <= 15) {
		int shift = i - 9;
		double lambda_1 = lambda(edge_lambdas[shift][0], p_glob);
		double lambda_2 = lambda(edge_lambdas[shift][1], p_glob);
		double ker_val = kernel(1, lambda_1 - lambda_2);
		result = lambda_1 * lambda_2 * ker_val;
	}
	// Функции третьего порядка, ассациированные с гранями
	else if (i <= 19) {
		int shift = i - 15;
		array<double, 3> face_lambda;
		result = 1;
		for(int f_i = 0; f_i < 3; f_i++) {
			face_lambda[i] = lambda(face_lambdas[shift][f_i], p_glob);
			result *= face_lambda[i];
		}
	}

	return result; 
}

vec3d tetelement::scalar_basis_grad_v(int i, double x, double y, double z) {
	point p_glob(x, y, z);
	vec3d result;



	// То, что заместо первого порядка
	if (i <= 3) {
		double lambda_i = lambda(i, p_glob);
		vec3d lambda_g = grad_lambda(i);
		result =  (2*lambda_i - 0.5) * lambda_g;
	}
	// Рёберные функции второго порядка
	else if(i <= 9) {
		int shift = i - 3;
		double lambda_1 = lambda(edge_lambdas[shift][0], p_glob);
		double lambda_2 = lambda(edge_lambdas[shift][1], p_glob);
		vec3d lambda_1_g = grad_lambda(edge_lambdas[shift][0]);
		vec3d lambda_2_g = grad_lambda(edge_lambdas[shift][1]);
		result = 2 * (lambda_1 * lambda_2_g + lambda_2 * lambda_1_g);
	}
	// Рёберные функции третьего порядка
	else if (i <= 15) {
		int shift = i - 9;
		double lambda_1 = lambda(edge_lambdas[shift][0], p_glob);
		double lambda_2 = lambda(edge_lambdas[shift][1], p_glob);
		double ker_val = kernel(1, lambda_1 - lambda_2);

		vec3d lambda_1_g = grad_lambda(edge_lambdas[shift][0]);
		vec3d lambda_2_g = grad_lambda(edge_lambdas[shift][1]);
		double d_ker = kernel_d(1, lambda_1 - lambda_2);

		result = lambda_1 * ker_val * lambda_2_g + lambda_2 * ker_val* lambda_1_g + lambda_1 * lambda_2 * d_ker * (lambda_1_g - lambda_2_g);
	}
	// Функции третьего порядка, ассациированные с гранями
	else if (i <= 19) {
		int shift = i - 15;
		array<double, 3> face_lambda;
		array<vec3d, 3> face_grads;
		for(int f_i = 0; f_i < 3; f_i++) {
			face_lambda[i] = lambda(face_lambdas[shift][f_i], p_glob);
			face_grads[i] = grad_lambda(face_lambdas[shift][f_i]);
		}

		result = face_lambda[0] * face_lambda[1] * face_grads[2] + face_lambda[0] * face_lambda[2] * face_grads[1] + face_lambda[1] * face_lambda[2] * face_grads[0];
	}

	return result; 
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

	// Определим барицентр тетраэдра
	for(int i = 0; i < 3; i++) {
		center[i] = 0;
		for(int j = 0; j < 4; j++)
			center[i] += node_array[j][i] / 4.0;
	}
	

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

	double Gauss_cord_gl[4][gauss_points_tet];



	//Перевод на текущий тетраэдр

	for(int i = 0; i < 4; i++) {
		for(int j = 0 ; j < gauss_points_tet; j++) {
			Gauss_cord_gl[i][j] = 0;
			for(int k = 0; k < 4; k++)
				Gauss_cord_gl[i][j] += D_matrix[i][k] * tet_integration::gauss_points_master[j][k];
		}
	}

	for(int i = 0; i < 4; i++)
		for(int j = 0; j < 3; j++)
			gauss_points[i][j] = Gauss_cord_gl[j][i];

	double w_t = 0;
	for(int i = 0; i < gauss_points_tet; i++) {
		gauss_weights[i] = tet_integration::gauss_weights[i] / 6;
		w_t += gauss_weights[i];
	}

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
				return mu * scalar_basis_v(i,x,y,z) * scalar_basis_v(j,x,y,z);
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
			return rp_func(x,y,z) * scalar_basis_v(i,x,y,z);
	});

	return b_loc;
}

double tetelement::L2_diff(func3d f, vector<double>& q_loc){

	int loc_n = q_loc.size();

	double res = integrate([&](double x, double y, double z){
		double u = 0;
		for (int j = 0; j < loc_n; j++) {
			double basis_v = scalar_basis_v(j, x, y, z);
			u += q_loc[j] * basis_v;
		}

		double f_v = f(x, y, z);
		return (u - f_v)*(u - f_v);
	});

	return res;
}

double tetelement::L2_diff(func3d f, vector<double>& q_loc, vector<double>& q_virtual){

	int loc_n = q_loc.size();
	int virt_n = q_virtual.size();

	double res = integrate([&](double x, double y, double z){
		double u = 0;
		for (int j = 0; j < loc_n; j++) {
			double basis_v = scalar_basis_v(j, x, y, z);
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