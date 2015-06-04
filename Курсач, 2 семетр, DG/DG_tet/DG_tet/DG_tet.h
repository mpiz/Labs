#include "BaseElementDG.h"

double solution(double x, double y, double z);

class DG_tet : public BaseElement<tetelement, trface> {
public:
	DG_tet();
	void calculate();
	void input_mesh(string file_name);

	double diff_L2(func3d f);
private:

	int get_order(vector<node>& nodes_s);
	double get_lambda(tetelement& el);

	map<tuple<int, int, int>, int> face_map;


	void make_face(tetelement* el, array<int, 3>& nodes_numbers);
	void make_faces();
	
};