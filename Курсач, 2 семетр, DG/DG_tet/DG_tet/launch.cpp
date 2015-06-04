#include "DG_tet.h"

int main() {
	DG_tet DG;

	DG.input_mesh("mesh_test.dat");
	DG.calculate();

	double error = DG.diff_L2(solution);

	return 0;
}