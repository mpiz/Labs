#include "VFEM_o2_t1_quad.h"
#include <tuple>



class VFEM_o2_t1_cube: public BaseElement<cubeelement> {
 public:
	VFEM_o2_t1_cube();

	double get_lambda(cubeelement& el);

	void input_mesh(string file_nodes, string file_elements);
	void input_bound(string file_elements);

	void calculate();

	void output_weights(string file_name, string file_bound);

	cubeelement* find_element(point pn);

 private:

	dof_type get_dof_number(vector<node>& el_nodes);
	void create_element(vector<int>& nodes_numbs, int el_i);
	void create_bound_element(vector<int>& nodes_numbs, int el_i);
	void prepare_boundry();

	enum class GEOM_TYPES{EDGE, FACE, VOLUME};

	dof_type get_geom_dof(GEOM_TYPES geom, int n_g, int num);

	int add_edge(int n1, int n2);
	int add_face(int n1, int n2, int n3, int n4);

	map<tuple<int, int>, int> edges; // Массив номеров рёбер
	map<tuple<int, int, int, int>, int> faces; // Массив номеров граней

	unordered_map<string, dof_type> dofs_map; // Список глобальных степеней свободы и соответвие их элементу

	dof_type cur_dof_count;

	int cur_edge_count, cur_face_count;

	VFEM_o2_t1_quad boundory_method;

	vector<dof_type> first_bound_dofs;
	set<dof_type> first_bound_dofs_set;
	vector<quadelement> bound_faces;
	int first_bound_n;
	int first_bound_elements_n;



};