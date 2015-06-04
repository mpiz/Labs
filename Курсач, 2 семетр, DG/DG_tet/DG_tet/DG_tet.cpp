#include "DG_tet.h"

double solution(double x, double y, double z) {
	return x*x*x;
}

double right_part(double x, double y, double z) {
	return -6*x;
}

DG_tet::DG_tet() {
	dofs_n = 1;
}

void DG_tet::input_mesh(string file_name) {
	BaseElement::input_mesh(file_name, 304);
	make_faces();
}

void DG_tet::make_face(tetelement* el, array<int, 3>& nodes_numbers) {
	trface* face;
	// Проверим - есть ли уже данная грань
	auto face_key = make_tuple(nodes_numbers[0], nodes_numbers[1], nodes_numbers[2]);

	auto face_iter = face_map.find(face_key);
	if (face_iter != face_map.end()) {
		face = &elements_faces[face_iter->second];
	}
	else {
		vector<node> nodes_s;
		for(int i = 0; i < 3; i++)
			nodes_s.push_back(local_nodes[nodes_numbers[i]]);
		int size = elements_faces.size();
		elements_faces.push_back(trface(nodes_s));
		face = &elements_faces[size];
		face_map[face_key] = size;

	}

	face->add_element(el);

}

void DG_tet::make_faces() {
	
	array<int, 3> nodes_numbers;
	for(int el_i = 0; el_i < elements_n; el_i++) {
		//Соберём все возможные грани с тетраэдра
		for(int i = 0; i <4; i++) {
			for(int j = 0; j < i; j++) {
				for(int k = 0; k < j; k++) {
					nodes_numbers[0] = elements[el_i].get_local_node(i).number;
					nodes_numbers[1] = elements[el_i].get_local_node(j).number;
					nodes_numbers[2] = elements[el_i].get_local_node(k).number;

					make_face(&elements[el_i], nodes_numbers);
				}
			}
		}
	}

	faces_n = elements_faces.size();
}


void DG_tet::calculate() {
	vector<func3d> rps;
	rps.push_back(right_part);
	rps.push_back(solution);

	generate_port();
	generate_matrix_with_out_bound(rps);

	solver.init(gi.data(), gj.data(), gu.data(), gl.data(), di.data(), rp[0], local_dof_n);
	solver.solve(solutions[0]);

}

double DG_tet::diff_L2(func3d f) {
	double diff = 0;
	double f_norm = 0;
	double* sol = solutions[0];
	for(int i = 0; i < elements_n; i++) {
		auto element_dofs = elements[i].get_dofs();
		vector<double> q_loc;
		for(auto iter = element_dofs.begin(); iter != element_dofs.end(); iter++)
			q_loc.push_back(sol[*iter]);

		diff += elements[i].L2_diff(f, q_loc);
		f_norm += elements[i].integrate([&](double x, double y, double z)->double { return f(x,y,z) * f(x,y,z);});
	}

	diff = sqrt(diff);
	f_norm = sqrt(f_norm);
	return diff;
}

double DG_tet::get_lambda(tetelement& el) {
	return 1.0;
}

int DG_tet::get_order(vector<node>& nodes_s) {
	return 3;
}