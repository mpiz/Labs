#include "DG_tet.h"

int main() {
	DG_tet DG;

	DG.input_mesh("mesh_3.dat");
	DG.calculate();

	double error = DG.diff_L2([&](double x, double y, double z)->double {
		return 1.0;
	});

	return 0;
}