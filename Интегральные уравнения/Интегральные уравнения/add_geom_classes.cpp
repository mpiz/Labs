#include "add_geom_classes.h"

plane::plane() {
	jacobian = 0;
}

plane::plane(const array<node, 3>& points) {
	plane_points = points;
	init_cords();
}

double plane::get_jacobian() const {
	return jacobian;
}

vec3d plane::get_base_vec(size_t i) const {
	return base_vec[i];
}

vec3d plane::get_tau(size_t i) const {
	return tau[i];
}

matrix(3) plane::get_matrix() const{
	return transition_matrix;
}

array<point,3> plane::get_tr_point() {
	return trpoint;
}

vec3d plane::get_normal() const {
	return normal_vector;
}

vec3d plane::get_normal_in_plane(vec3d a) const {
	matrix(3) turn_matrix;
	turn_matrix[0][0] = 0.0; turn_matrix[0][1] = -1.0; turn_matrix[0][2] = 0.0;
	turn_matrix[1][0] = 1.0; turn_matrix[1][1] = 0.0; turn_matrix[1][2] = 0.0;
	turn_matrix[2][0] = 0.0; turn_matrix[2][1] = 0.0; turn_matrix[2][2] = 0.0;

	vec3d a_loc = to_local_cord(a);
	vec3d a_loc_n = turn_matrix * a_loc;
	vec3d a_n = to_global_cord(a_loc_n);

	return a_n;

	
}

point plane::to_local_cord(point p_glob) {
	point p_shift = p_glob;

	//�����

	for(int i = 0; i < 3; i++)
		p_shift[i] -= plane_points[0][i];

	point p_loc;

	//�������
	for(int i = 0; i < 3; i++) {
		p_loc[i] = 0;
		for(int j = 0; j < 3; j++)
			p_loc[i] += transition_matrix[i][j] * p_shift[j];
	}

	return p_loc;
}

vec3d plane::to_local_cord(vec3d v_glob) const {

	vec3d v_loc = transition_matrix * v_glob;
	return v_loc;
}


point plane::to_global_cord(point p_loc) {
	point p_turn;

	//�������
	for(int i = 0; i < 3; i++) {
		p_turn[i] = 0;
		for(int j = 0; j < 3; j++)
			p_turn[i] += transition_matrix[j][i] * p_loc[j];
	}

	point p_glob; 

	//�����
	for(int i = 0; i < 3; i++)
		p_glob[i] = p_turn[i] + plane_points[0][i];

	return p_glob;

}

vec3d plane::to_global_cord(vec3d v_loc) const {
	vec3d v_glob;

	//�������
	for(int i = 0; i < 3; i++) {
		v_glob[i] = 0;
		for(int j = 0; j < 3; j++)
			v_glob[i] += transition_matrix[j][i] * v_loc[j];
	}


	return v_glob;

}

void plane::init_cords() {
		//���������� ��������� ������� ���������
	base_vec[0] = vec3d(plane_points[0], plane_points[1]); 
	base_vec[1] = vec3d(plane_points[0], plane_points[2]); 
	vec3d e2 = base_vec[1] - ((base_vec[0]*base_vec[1]) / (base_vec[0]*base_vec[0])) * base_vec[0];
	vec3d e3 = base_vec[0].cross(e2);
	normal_vector = e3;

	tau[0] = base_vec[0];
	tau[1] = base_vec[1];
	tau[2] = vec3d(plane_points[1], plane_points[2]);

	h[0] = base_vec[0].norm();
	h[1] = base_vec[1].norm();

	trpoint[0] = point(0, 0, 0);
	trpoint[1] = point(base_vec[0].norm(), 0, 0);

	vec3d e1 = base_vec[0] / base_vec[0].norm();
	e2 = e2 / e2.norm();
	e3 = e3 / e3.norm();

	for(int i = 0; i < 3; i++) {
		transition_matrix[0][i] = e1[i].real();
		transition_matrix[1][i] = e2[i].real();
		transition_matrix[2][i] = e3[i].real();
	}

	trpoint[2] = (transition_matrix * base_vec[1]).to_point();
	jacobian = (base_vec[0].cross(base_vec[1])).norm();
}