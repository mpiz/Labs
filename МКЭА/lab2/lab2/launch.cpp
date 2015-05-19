#include "VFEM_o2_t1_cube.h"
#include "gen_grid.h"
#include <stdio.h>

int main() {

	generate_3unreg_grid();
	generate_ku_faces(2, 2, 2);

	VFEM_o2_t1_cube fem;
	fem.input_mesh("xyz.txt", "nvtr_v.txt");
	fem.input_bound("ku_faces.txt");

	fem.calculate();

	fem.output_weights("solution.txt", "bound.txt");

	return 0;
}