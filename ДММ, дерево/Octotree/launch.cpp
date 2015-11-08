#include "VFEM_E.h"
#include <stdlib.h>
#include <stdio.h>


int main() {


	VFEM_E octtree;
	octtree.input_mesh("web.msh");

	int samples_n = 1e4;
	double x, y, z;
	double tree_time = 0, linear_time = 0;
	double tt_iter, lt_iter;

	ofstream outp_info("iters_info.txt");

	srand(GetTickCount());

	for (int i = 0; i < samples_n; i++) {

		x = rand()%100 * 1.0 / 100;
		y = rand()%100 * 1.0 / 100;
		z = rand()%100 * 1.0 / 100;

		tt_iter = octtree.function_in_point_tree(x, y, z);
		lt_iter = octtree.function_in_point_linear(x, y, z);

		tree_time += tt_iter;
		linear_time += lt_iter;

		outp_info << i << "\t" << x << "\t" << y << "\t" << z << "\t" << tt_iter << "\t" << lt_iter << endl;
	}

	printf("Avg time, tree method: %.3e, Full time: %.3e\n", tree_time / samples_n, tree_time);
	printf("Avg time, linear method: %.3e, Full time: %.3e\n", linear_time / samples_n, linear_time);
	outp_info.close();

	system("pause");

	return 0;
}