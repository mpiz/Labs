#include "BaseElement.h"

class VFEM_o2_t1_quad : public BaseElement<quadelement> {
 public:
	 void set_elements(vector<quadelement>& s_faces, set<dof_type>& bound_dofs);
	 void calculate();
 private:


};