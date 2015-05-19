#include "VFEM_o2_t1_quad.h"

double VFEM_o2_t1_quad::get_dof_weight(dof_type glob_dof_i) {
	auto loc_dof = glob_to_loc[glob_dof_i];

	return solutions[0][loc_dof];
}

void VFEM_o2_t1_quad::set_elements(vector<quadelement>& s_faces, set<dof_type>& bound_dofs) {
	// Инициализируем степени свободы
	dofs_n = 1;
	local_dof_n = 0;
	for(auto iter = bound_dofs.begin(); iter != bound_dofs.end(); iter++, local_dof_n++) {
		glob_to_loc[*iter] = local_dof_n;
		loc_to_glob[local_dof_n] = *iter;
	}

	//Добавим элементы
	elements_n = s_faces.size();
	elements.resize(elements_n);
	for(int el_i = 0; el_i < elements_n; el_i++) {
		elements[el_i] = s_faces[el_i];
		elements[el_i].renumerate_dofs(glob_to_loc);
	}

}

void VFEM_o2_t1_quad::calculate() {
	generate_port();
	vector<vfunc3d> eq_rp;
	eq_rp.push_back(
		[](double x, double y, double z)->vec3d {
		return vec3d(1.0, 0, 0);	
		}
	);
	generate_matrix_with_out_bound(eq_rp);

	int basis_i = 0;

#ifdef DEBUGOUTP
	auto rp_v = rp[basis_i];
	auto sol = solutions[basis_i];
#endif

	solver.init(&gi.front(), &gj.front(), &di.front(), &gg.front(), local_dof_n);
	solver.solve(rp[basis_i], solutions[basis_i]);


}

void VFEM_o2_t1_quad::output_weights(string file_name) {
	ofstream outp_file(file_name.c_str());

	for(int i = 0; i < local_dof_n; i++)
		outp_file << i << "\t" << loc_to_glob[i] << "\t" <<solutions[0][i] << endl;

	outp_file.close();

}

quadelement* VFEM_o2_t1_quad::find_element(point pn) {
	return nullptr;
}

double VFEM_o2_t1_quad::get_lambda(quadelement& el){
	return 1.0;
}

