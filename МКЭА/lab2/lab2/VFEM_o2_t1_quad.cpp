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
		return vec3d(x, 0, 0);	
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

double VFEM_o2_t1_quad::get_lambda(quadelement& el){
	return 1.0;
}

