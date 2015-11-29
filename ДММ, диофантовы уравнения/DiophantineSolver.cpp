#include "DiophantineSolver.h"

void DiophantineSolver::input(string file_name) {
	ifstream inp_file(file_name);

	inp_file >> m >> n;

	full_matrix_n = n + 1;
	full_matrix_m = m + n;
	matrix_A.resize(full_matrix_m);
	for(int i = 0; i < m; i++) {
		matrix_A[i].resize(full_matrix_n);
		for(int j = 0; j < full_matrix_n; j++)
			inp_file >> matrix_A[i][j];
		matrix_A[i][n] *= -1;
	}

	// Допишем единичную
	for(int i = m; i < full_matrix_m; i++) {
		matrix_A[i].resize(full_matrix_n);
		for(int j = 0; j < full_matrix_n; j++)
			matrix_A[i][j] = (i-m) == j ? 1 : 0;
	}

	inp_file.close();

}

void DiophantineSolver::output(string file_name) {
	ofstream outp_file(file_name);

	if (solution_result) {

		outp_file << s << endl;
		for(int i = 0; i < n; i++) {
			// В каждой строке s+1 число: частное решение и s коэфф-тов при свободных переменных
			outp_file << matrix_A[m + i][n];
			if(s > 0) {
				outp_file << " ";
				for(int j = 0; j < s - 1; j++)
					outp_file << matrix_A[m + i][n - s + j] << " ";
				outp_file << matrix_A[m + i][n - 1] << endl;
			}
			else
				outp_file << endl;
		}
	}
	else
		outp_file << "Error!" << endl;

	outp_file.close();

}

int DiophantineSolver::find_min_in_row(int i) {
	int min_ind = -1;
	int min_val = numeric_limits<int>::max();
	for(int j = i; j < n; j++) {
		int cur_val = abs(matrix_A[i][j]);
		if(min_val > cur_val && cur_val != 0) {
			min_ind = j;
			min_val = abs(matrix_A[i][j]);
		}
	}
	return min_ind;
}

void DiophantineSolver::transform_column(int i, int min_j, int val) {
	for(int j = i; j < full_matrix_n; j++) {
		if (j != min_j && matrix_A[i][j] != 0) {
			int k = matrix_A[i][j] / val;
			for(int i1 = i; i1 < full_matrix_m; i1++) 
				matrix_A[i1][j] -= k * matrix_A[i1][min_j];
		}

	}

	if (matrix_A[i][i] == 0) {
		for(int i1 = i; i1 < full_matrix_m; i1++)
			matrix_A[i1][i] += matrix_A[i1][min_j];

	}

}

bool DiophantineSolver::check_status(int i) {
	bool st = true;
	for(int j = i + 1; j < n; j++) {
		if (matrix_A[i][j] != 0) {
			st = false;
			exit;
		}
	}

	if (st) {

		if ( matrix_A[i][n] % matrix_A[i][i] != 0)
			solution_result = false;			
		st = (matrix_A[i][n] == 0);
	}

	return st;
}

bool DiophantineSolver::solve() {

	solution_result = true;

	for(int i = 0; i < m && solution_result; i++) {
		bool row_in_procces = true;
		while (row_in_procces && solution_result) {
			// Найдём миниальный по модулю, отличный от нуля элемент строки
			int min_j = find_min_in_row(i);

			// Если строка нулевая, то выкинем её из СЛАУ и перейдём на следующую или выдадаим ошибку
			if (min_j == -1) {
				if (matrix_A[i][n] == 0) {
					matrix_A.erase(matrix_A.begin() + i);
					m--;
					full_matrix_m--;
					i--;
				}
				else
					solution_result = false;
				break;
			}

			transform_column(i, min_j, matrix_A[i][min_j]);
			row_in_procces = !check_status(i);

		}

	}
	if (solution_result) {

		s = n - m;
		if (s < 0) {
			solution_result = false;
		}

	}
	return solution_result;
}

