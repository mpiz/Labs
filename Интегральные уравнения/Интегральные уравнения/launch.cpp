#include "IntegralEq.h"

int main () {
	IntegralEq inteq;

	inteq.set_env(10, 0.1, 1, point(-100, -100, -10), point(100, 100, -2));
	inteq.input_mesh("Group_1.dat");
	cout << "Reading mesh: completed\n";
	inteq.calculate();

	cout << "Output result\n";
	inteq.output("res.txt");

	cout << "Output E0\n";
	inteq.outputE0("res_E0.txt");

	cout << "Output surface\n";
	inteq.outputE_surface("res_surf.txt");
	return 0;
}