#include "integration_points.h"

namespace integration_points {
		
		const double QUAD_7_a = sqrt((144.0-3.0*sqrt(583.0))/287.0);
		const double QUAD_7_b = sqrt((144.0+3.0*sqrt(583.0))/287.0);
		const double QUAD_7_c = sqrt(6.0 / 7.0);

		const double QUAD_7_wa = 307.0 / 810.0 + 923.0 / (270.0 * sqrt(583.0));
		const double QUAD_7_wb = 307.0 / 810.0 - 923.0 / (270.0 * sqrt(583.0));
		const double QUAD_7_wc = 98.0 / 405.0;

		const int QUAD_7_n = 12;
		 
		integ_point QUAD_7_p[QUAD_7_n] = {
			integ_point(-QUAD_7_c, 0),
			integ_point(QUAD_7_c, 0),
			integ_point(0, -QUAD_7_c),
			integ_point(0, QUAD_7_c),
			integ_point(-QUAD_7_a, -QUAD_7_a),
			integ_point(QUAD_7_a, -QUAD_7_a),
			integ_point(-QUAD_7_a, QUAD_7_a),
			integ_point(QUAD_7_a, QUAD_7_a),
			integ_point(-QUAD_7_b, -QUAD_7_b),
			integ_point(QUAD_7_b, -QUAD_7_b),
			integ_point(-QUAD_7_b, QUAD_7_b),
			integ_point(QUAD_7_b, QUAD_7_b)
		};

		double QUAD_7_w[QUAD_7_n] = {
			QUAD_7_wc,
			QUAD_7_wc,
			QUAD_7_wc,
			QUAD_7_wc,
			QUAD_7_wa,
			QUAD_7_wa,
			QUAD_7_wa,
			QUAD_7_wa,
			QUAD_7_wb,
			QUAD_7_wb,
			QUAD_7_wb,
			QUAD_7_wb
		};


		const double CUBE_7_a = sqrt(6.0 / 7.0);
		const double CUBE_7_b = sqrt( (960.0 - 33.0*sqrt(238.0)) / 2726.0);
		const double CUBE_7_c = sqrt( (960.0 + 33.0*sqrt(238.0)) / 2726.0);

		const double CUBE_7_w_1 = 1078.0 / 3645.0;
		const double CUBE_7_w_2 = 343.0 / 3645.0;
		const double CUBE_7_w_3 = 43.0/135.0 + 829.0 * sqrt(238.0) / 136323.0;
		const double CUBE_7_w_4 = 43.0/135.0 - 829.0 * sqrt(238.0) / 136323.0;

		const int CUBE_7_n = 34;

		integ_point CUBE_7_p[CUBE_7_n] = {
			integ_point(CUBE_7_a, 0, 0),
			integ_point(-CUBE_7_a, 0, 0),
			integ_point(0, CUBE_7_a, 0),
			integ_point(0, -CUBE_7_a, 0),
			integ_point(0, 0, CUBE_7_a),
			integ_point(0, 0, -CUBE_7_a),

			integ_point(CUBE_7_a, CUBE_7_a, 0),
			integ_point(CUBE_7_a, -CUBE_7_a, 0),
			integ_point(-CUBE_7_a, CUBE_7_a, 0),
			integ_point(-CUBE_7_a, -CUBE_7_a, 0),
			integ_point(CUBE_7_a, 0, CUBE_7_a),
			integ_point(CUBE_7_a, 0, -CUBE_7_a),
			integ_point(-CUBE_7_a, 0, CUBE_7_a),
			integ_point(-CUBE_7_a, 0, -CUBE_7_a),
			integ_point(0, CUBE_7_a, CUBE_7_a),
			integ_point(0, CUBE_7_a, -CUBE_7_a),
			integ_point(0, -CUBE_7_a, CUBE_7_a),
			integ_point(0, -CUBE_7_a, -CUBE_7_a),

			integ_point(CUBE_7_b, CUBE_7_b, CUBE_7_b),
			integ_point(CUBE_7_b, CUBE_7_b, -CUBE_7_b),
			integ_point(CUBE_7_b, -CUBE_7_b, CUBE_7_b),
			integ_point(CUBE_7_b, -CUBE_7_b, -CUBE_7_b),
			integ_point(-CUBE_7_b, CUBE_7_b, CUBE_7_b),
			integ_point(-CUBE_7_b, CUBE_7_b, -CUBE_7_b),
			integ_point(-CUBE_7_b, -CUBE_7_b, CUBE_7_b),
			integ_point(-CUBE_7_b, -CUBE_7_b, -CUBE_7_b),

			integ_point(CUBE_7_c, CUBE_7_c, CUBE_7_c),
			integ_point(CUBE_7_c, CUBE_7_c, -CUBE_7_c),
			integ_point(CUBE_7_c, -CUBE_7_c, CUBE_7_c),
			integ_point(CUBE_7_c, -CUBE_7_c, -CUBE_7_c),
			integ_point(-CUBE_7_c, CUBE_7_c, CUBE_7_c),
			integ_point(-CUBE_7_c, CUBE_7_c, -CUBE_7_c),
			integ_point(-CUBE_7_c, -CUBE_7_c, CUBE_7_c),
			integ_point(-CUBE_7_c, -CUBE_7_c, -CUBE_7_c)

		};
		double CUBE_7_w[CUBE_7_n] = {
			CUBE_7_w_1, CUBE_7_w_1, CUBE_7_w_1, CUBE_7_w_1,
			CUBE_7_w_1, CUBE_7_w_1,
			CUBE_7_w_2, CUBE_7_w_2, CUBE_7_w_2, CUBE_7_w_2,
			CUBE_7_w_2, CUBE_7_w_2, CUBE_7_w_2, CUBE_7_w_2,
			CUBE_7_w_2, CUBE_7_w_2, CUBE_7_w_2, CUBE_7_w_2,
			CUBE_7_w_3, CUBE_7_w_3, CUBE_7_w_3, CUBE_7_w_3,
			CUBE_7_w_3, CUBE_7_w_3, CUBE_7_w_3, CUBE_7_w_3,
			CUBE_7_w_4, CUBE_7_w_4, CUBE_7_w_4, CUBE_7_w_4,
			CUBE_7_w_4, CUBE_7_w_4, CUBE_7_w_4, CUBE_7_w_4
		};

};


void get_QUAD(int_type order, integ_point* &points, double* &weights, int_type& n) {

	switch(order) {
		case 7:
			points = integration_points::QUAD_7_p;
			weights = integration_points::QUAD_7_w;
			n = integration_points::QUAD_7_n;
			break;

	default:
		throw ERRS::ORDER;
	}
}


void get_CUBE(int_type order, integ_point* &points, double* &weights, int_type& n) {

	switch(order) {
		case 7:
			points = integration_points::CUBE_7_p;
			weights = integration_points::CUBE_7_w;
			n = integration_points::CUBE_7_n;
			break;

	default:
		throw ERRS::ORDER;
	}
}



void get_quadrature_rule(ELEM_TYPE element_type, int_type order, integ_point* &points, double* &weights, int_type& n) {
	switch(element_type) {
		case ELEM_TYPE::QUAD :
			get_QUAD(order, points, weights, n);
			break;
		case ELEM_TYPE::CUBE :
			get_CUBE(order, points, weights, n);
			break;
		default:
			throw ERRS::TYPE;

	}
}

