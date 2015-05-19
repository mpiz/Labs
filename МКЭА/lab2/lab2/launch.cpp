#include "VFEM_o2_t1_cube.h"
#include "gen_grid.h"
#include <stdio.h>

int main() {

	generate_3unreg_grid();
	generate_ku_faces(1, 1, 1);

	VFEM_o2_t1_cube fem;
	fem.input_mesh("xyz.txt", "nvtr_v.txt");
	fem.input_bound("ku_faces.txt");

	fem.calculate();

	fem.output_weights("solution.txt");

	/*for(int i = 0; i < 39; i++)
		dofs.push_back(i);

	cubeelement cube(points, dofs);

	double res = cube.integrate([](double x, double y, double z)->double { return pow(x, 6);});


	dyn_matrix M = cube.get_local_matrix(1.0, 0.0);

	FILE* M_f = fopen("M.txt", "w");
	int m_s = M.size();
	for(int i = 0; i < m_s; i++) {
		for(int j = 0; j < m_s; j++) {
			fprintf(M_f, "%lf", M[i][j]);
			if (j == m_s - 1) 
				fprintf(M_f, "%\n", M[i][j]);
			else
				fprintf(M_f, "%\t", M[i][j]);
		}

	}

	fclose(M_f);*/

	return 0;
}