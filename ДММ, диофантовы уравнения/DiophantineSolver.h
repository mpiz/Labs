#include <vector>
#include <string>
#include <fstream>
#include <utility>
#include <limits>

using namespace std;

class DiophantineSolver {
public:

	void input(string file_name);
	bool solve();
	void output(string file_name);

private:
	int n, m, s;
	vector< vector<int> > matrix_A;

	int full_matrix_n, full_matrix_m;

	int find_min_in_row(int i);
	void transform_column(int i, int min_j, int val);
	bool check_status(int i);

	bool solution_result;

	

};