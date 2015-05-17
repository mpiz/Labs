#pragma once

#include <cmath>

/*
		���� ��� ������ ������ ������������ ������ ��� ������ ���������.
	������� get_quadtature_rule ���������� ����� �������������� ��� �������� �������� � �������. ���� ����� ��� - ������:
		element_type (in) - ��� ��������;
		order (in) - ������� ��������������;
		points (out) - ��������� ����� ��������������;
		weights (out) - ���� ��������������;
		n (out) - ���������� ����� ��������������.

	���� ��������� �������������� � ��������� �������:
		QUAD - ������� [-1, 1] ������-������� ��� ����������������.
			7
		CUBE - ��� [-1, 1] ������-������� ��� ��������������
			7
*/

typedef int int_type;

struct integ_point {
	integ_point() {

	}

	integ_point(double x) {
		xi = x;
	}

	integ_point(double x, double y) {
		xi = x; eta = y;
	}
	integ_point(double x, double y, double z) {
		xi = x; eta = y; zeta = z;
	}
	double xi;
	double eta;
	double zeta;
};

enum class ERRS{TYPE, ORDER};
enum class ELEM_TYPE {QUAD, CUBE};


void get_quadrature_rule(ELEM_TYPE element_type, int_type order, integ_point* &points, double* &weights, int_type& n); 
