#include "VFEM_E.h"

#include "stdio.h"

const int triangle_el = 2;
const int tethedron_el = 4;

// == Основные функции ==

void VFEM_E::input_mesh(string inp_file) {
	ifstream inp_f(inp_file.c_str());

	//Считывание узлов
	inp_f >> nodes_n;

	node read_n;

	for(int i = 0; i < nodes_n; i++) {
		inp_f >> read_n;
		read_n.number--;
		nodes.push_back(read_n);
	}

	//Пролистываем
	string read_str;
	while(read_str != "$Elements")
		inp_f >> read_str;

	//Считывание граней и элементов

	int total_els, tmp_int, el_type;
	inp_f >> total_els;
	int ph_area;
	int params_n;
	array<int, 3> trnode;
	array<int, 3> tredge;
	array<int, 4> tetnode;
	array<int, 6> tetedge;

	elements.reserve(total_els);
	
	for(int iter = 0; iter < total_els; iter++) {
		inp_f >> tmp_int >> el_type;
		inp_f >> params_n;
		inp_f >> ph_area;

		if((iter * 100 / total_els) % 5 == 0) {
			cout << iter << "\\" << total_els << "\r";
			cout.flush();
		}


		switch(el_type) {
		case triangle_el: {
			for(int i = 0; i < params_n-1; i++);
				inp_f >> tmp_int;
	
			for(int i = 0; i < 3; i++) {
				inp_f >> trnode[i];
				trnode[i]--;
			}
		} break;

		case tethedron_el: {
			for(int i = 0; i < params_n-1; i++);
				inp_f >> tmp_int;
	
			for(int i = 0; i < 4; i++) {
				inp_f >> tetnode[i];
				tetnode[i]--;
			}

			sort(tetnode.begin(), tetnode.end());

			elements.push_back(tetelement(nodes[tetnode[0]], nodes[tetnode[1]], nodes[tetnode[2]], nodes[tetnode[3]]));
			int add_num = elements.size() - 1;
			elements[add_num].set_ph_area(ph_area);

			} break;

		}; 

	}

	elements_n = elements.size();

	//Построение дерева поиска

	ifstream gabs("gabarits.txt");
	double gx0, gx1, gy0, gy1, gz0, gz1;

	gabs >> gx0 >> gx1 >> gy0 >> gy1 >> gz0 >> gz1;

	gabs.close();

	cout << "Search tree construct begin" << endl;

	search_tree.construct_tree(gx0, gx1, gy0, gy1, gz0, gz1, elements);
	cout << "Search tree construct complete" << endl;
}


vec3d VFEM_E::function_in_point(double x, double y, double z) {

	tetelement* point_el = search_tree.find_point(x, y, z);
	if(point_el != NULL) {
		int loc_dof[12];
			for(int i = 0; i < 6; i++) {
				loc_dof[i] = (*point_el)[i];
				loc_dof[i+6] = (*point_el)[i] + edges_n;
			}
			vec3d val(0,0,0);

			for(int i = 0; i < 12; i++) {
				vec3d basis_v = point_el->basis_v(i, x, y, z);
				val = val + solution[loc_dof[i]] * basis_v;
			}

			return val;

	}
	return vec3d(0,0,0);

}

