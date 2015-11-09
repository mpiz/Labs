#include "elements_classes.h"


// ======== Треугольники ========

trelement::trelement() {
}

trelement::trelement(node n1, node n2, node n3) {
	node_array[0] = n1;
	node_array[1] = n2;
	node_array[2] = n3;

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

void trelement::init_cords() {

	basis[0] = &trelement::w1;
	basis[1] = &trelement::w2;
	basis[2] = &trelement::w3;
	basis[3] = &trelement::w4;
	basis[4] = &trelement::w5;
	basis[5] = &trelement::w6;

	kernel_basis[0] = &trelement::gradphi1;
	kernel_basis[1] = &trelement::gradphi2;
	kernel_basis[2] = &trelement::gradphi3;
	kernel_basis[3] = &trelement::w4;
	kernel_basis[4] = &trelement::w5;
	kernel_basis[5] = &trelement::w6;

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

	for(int i = 0; i < 3; i++) {
		transition_matrix[0][i] = e1[i].real();
		transition_matrix[1][i] = e2[i].real();
		transition_matrix[2][i] = e3[i].real();
	}

	trpoint[2] = (transition_matrix * g2).to_point();

	generate_L_cords();

	//Нахождене точек интегрирования по Гауссу, в локальной системе координат

	jacobian = (g1.cross(g2)).norm();

	gauss_points[0] = point(h1 / 2, 0, 0);
	gauss_points[1] = point(trpoint[2][0] / 2.0, trpoint[2][1] / 2.0, trpoint[2][2] / 2.0);
	gauss_points[2] = (transition_matrix * (0.5*g1 + 0.5*g2)).to_point();

	gauss_weights[0] = gauss_weights[1] = gauss_weights[2] = 1.0 / 6.0;

	for(int i = 0; i < gauss_points_tr; i++)
		gauss_points_global[i] = to_global_cord(gauss_points[i]);

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

vec3d trelement::basis_v(int i, double x, double y, double z) {
	point p_glob(x,y,z);
	return (this->*basis[i])(p_glob);

}

vec3d trelement::get_tau(int i) {
	return tau[i];
}

dcomplex trelement::integrate(func3d integ_func) {

	dcomplex res = 0;

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

matrix(6) trelement::get_M_matrix() {
	matrix(6) M;

	for(int i = 0; i < 6; i++)
		for(int j = 0; j <= i; j++) {
			M[i][j] = integrate([&](double x, double y, double z)->dcomplex {
				point p_glob(x,y,z);
				return (this->*basis[i])(p_glob) * (this->*basis[j])(p_glob);
			}).real();
			M[j][i] = M[i][j];
		}

	return M;
}

matrix(6) trelement::get_kernel_M_matrix() {
	matrix(6) M;

	for(int i = 0; i < 6; i++)
		for(int j = 0; j <= i; j++) {
			M[i][j] = integrate([&](double x, double y, double z)->dcomplex {
				point p_glob(x,y,z);
				return (this->*kernel_basis[i])(p_glob) * (this->*kernel_basis[j])(p_glob);
			}).real();
			M[j][i] = M[i][j];
		}

	return M;

}

array<dcomplex, 6> trelement::get_right_part(vfunc3d rp_func) {
	array<dcomplex, 6> b;

	for(int i = 0; i < 6; i++)
		b[i] = integrate([&](double x, double y, double z)->dcomplex {
			point p_glob(x,y,z);
			vec3d rp_v = rp_func(x,y,z);
			vec3d basis_v = (this->*basis[i])(p_glob);
			return  rp_v * basis_v;
	});

	return b;
}

array<dcomplex, 6> trelement::get_kernel_right_part(vfunc3d rp_func) {
	array<dcomplex, 6> b;

	for(int i = 0; i < 6; i++)
		b[i] = integrate([&](double x, double y, double z)->dcomplex {
			point p_glob(x,y,z);
			vec3d rp_v = rp_func(x,y,z);
			vec3d basis_v = (this->*kernel_basis[i])(p_glob);
			return  rp_v * basis_v;
	});

	return b;
}

vec3d trelement::grad_lambda(int i) {
	vec3d grad;
	for(int j = 0; j < 2; j++) {
		grad[j] = L_cord_matrix[i][j];
	}

	return to_global_cord(grad);
}

// == Базисные функции ==
vec3d trelement::w1(point p_glob) {

	point p_loc = to_local_cord(p_glob);
	double l_coef[2];

	l_coef[0] = lambda(0, p_loc);
	l_coef[1] = lambda(1, p_loc);

	vec3d grads_x[2];
	grads_x[0] = grad_lambda(1);
	grads_x[1] = grad_lambda(0);

	return l_coef[0]*grads_x[0] - l_coef[1]*grads_x[1];
}

vec3d trelement::w2(point p_glob) {
	double l_coef[2];
	point p_loc = to_local_cord(p_glob);

	l_coef[0] = lambda(0, p_loc);
	l_coef[1] = lambda(2, p_loc);

	vec3d grads_x[2];
	grads_x[0] = grad_lambda(2);
	grads_x[1] = grad_lambda(0);

	return l_coef[0]*grads_x[0] - l_coef[1]*grads_x[1];
}

vec3d trelement::w3(point p_glob) {
	double l_coef[2];
	point p_loc = to_local_cord(p_glob);

	l_coef[0] = lambda(1, p_loc);
	l_coef[1] = lambda(2, p_loc);

	vec3d grads_x[2];
	grads_x[0] = grad_lambda(2);
	grads_x[1] = grad_lambda(1);

	return l_coef[0]*grads_x[0] - l_coef[1]*grads_x[1];
}

vec3d trelement::w4(point p_glob) {
	double l_coef[2];
	point p_loc = to_local_cord(p_glob);

	l_coef[0] = lambda(0, p_loc);
	l_coef[1] = lambda(1, p_loc);

	vec3d grads_x[2];
	grads_x[0] = grad_lambda(1);
	grads_x[1] = grad_lambda(0);

	return l_coef[0]*grads_x[0] + l_coef[1]*grads_x[1];
}

vec3d trelement::w5(point p_glob) {
	double l_coef[2];
	point p_loc = to_local_cord(p_glob);

	l_coef[0] = lambda(0, p_loc);
	l_coef[1] = lambda(2, p_loc);

	vec3d grads_x[2];
	grads_x[0] = grad_lambda(2);
	grads_x[1] = grad_lambda(0);

	return l_coef[0]*grads_x[0] + l_coef[1]*grads_x[1];
}

vec3d trelement::w6(point p_glob) {
	double l_coef[2];
	point p_loc = to_local_cord(p_glob);

	l_coef[0] = lambda(1, p_loc);
	l_coef[1] = lambda(2, p_loc);

	vec3d grads_x[2];
	grads_x[0] = grad_lambda(2);
	grads_x[1] = grad_lambda(1);

	return l_coef[0]*grads_x[0] + l_coef[1]*grads_x[1];
}

vec3d trelement::gradphi1(point p_loc) {
	return grad_lambda(0);
}

vec3d trelement::gradphi2(point p_loc) {
	return grad_lambda(1);
}

vec3d trelement::gradphi3(point p_loc) {
	return grad_lambda(2);
}



// ======== Тетраэдры ========

tetelement::tetelement() {
}

tetelement::tetelement(node n1, node n2, node n3, node n4) {
	node_array[0] = n1;
	node_array[1] = n2;
	node_array[2] = n3;
	node_array[3] = n4;

	init_coords();
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

bool tetelement::in_element(double x, double y, double z) {
	point p_glob(x,y,z);
	for(int i = 0; i < 4; i++) {
		double L = lambda(i, p_glob);
		if(L > 1 || L < 0)
			return false;
	}

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