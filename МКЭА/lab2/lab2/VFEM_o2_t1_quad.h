#include "BaseElement.h"

class VFEM_o2_t1_quad : public BaseElement<quadelement> {
 public:
	double get_lambda(quadelement& el);
	 void set_elements(vector<quadelement>& s_faces, set<dof_type>& bound_dofs);
	 void calculate();
	 double get_dof_weight(dof_type glob_dof_i);
	 void output_weights(string file_name);

	 quadelement* find_element(point pn);

};