#include "BaseElementDG.h"


class DG_tet : public BaseElement<tetelement, trface> {
public:
	void calculate();
	void input_mesh(string file_name);
private:

};