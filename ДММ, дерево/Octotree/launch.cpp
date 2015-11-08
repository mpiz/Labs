#include "VFEM_E.h"
#include <stdlib.h>


int main() {


	VFEM_E octtree;
	octtree.input_mesh("web.msh");

	double x, y, z;

	x = rand()%100 * 0.1 / 100 - 0.05;
	y = rand()%100 * 0.04 / 100 - 0.2;
	z = rand()%100 * 0.02 / 100 - 0.01;

	cout << "Searching for:\nx=" << x << "\ny=" << y << "\nz=" << z << endl;
	cout << "Using tree: ";

	double tree_time = octtree.function_in_point_tree(x, y, z);
	cout << tree_time << "\nUsing linear method: ";
	double linear_time = octtree.function_in_point_linear(x, y, z);
	cout << linear_time;

	system("pause");

	return 0;
}