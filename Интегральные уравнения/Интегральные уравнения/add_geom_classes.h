#include "geom_classes.h"


class plane {
 public:
	 plane();
	 plane(const array<node, 3>& points);

	point to_local_cord(point p_glob);
	vec3d to_local_cord(vec3d v_glob) const;

	point to_global_cord(point p_loc);
	vec3d to_global_cord(vec3d v_loc) const;


	double get_jacobian() const;
	vec3d get_base_vec(size_t i) const;
	vec3d get_tau(size_t i) const;
	vec3d get_normal() const;
	
	vec3d get_normal_in_plane(vec3d a) const;


	matrix(3) get_matrix() const;
	array<point,3> get_tr_point();
 private:

	 void init_cords(); //���������� ��������� ���������
	 array<node, 3> plane_points; // ��� �����, ������������ ���������

	array<point,3> trpoint; //��������� ���������� ����� ������������
	matrix(3) transition_matrix; //������� �������� � ��������� ����������

	array<vec3d, 4> tau;

	vec3d normal_vector;

	double jacobian;
	double h[2];
	array<vec3d, 2> base_vec;


};