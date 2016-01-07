#include "IntegralEq.h"

int main () {
	IntegralEq inteq;

	inteq.set_env(5, 0.1, 1, point(-2, -2, -5), point(2, 2, -2));
	inteq.input_mesh("Group_2.dat");
	cout << "Reading mesh: completed\n";
	inteq.calculate();

	//cout << "Output result\n";
	//inteq.output("res.txt");

	//cout << "Output E0\n";
	//inteq.outputE0("res_E0.txt");

	cout << "Output surface\n";
	inteq.outputE_surface_tp("res_surf0.txt");
	return 0;
}