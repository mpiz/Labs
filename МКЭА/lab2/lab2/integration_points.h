#pragma once

#include <cmath>

/*
		Файл имён хранит наборы квадратруных формул для мастер элементов.
	Функция get_quadtature_rule возвращает точки интегрирования для заданого элемента и порядка. Если таких нет - ошибка:
		element_type (in) - тип элемента;
		order (in) - порядок интегрирования;
		points (out) - множество точек интегрирования;
		weights (out) - веса интегрирования;
		n (out) - количество точек интегрирования.

	Типы элементов интегрирования и доступные порядки:
		QUAD - квадрат [-1, 1] мастер-элемент для четырёхугольников.
			7
		CUBE - куб [-1, 1] мастер-элемент для шестигранников
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
