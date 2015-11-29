#include "DiophantineSolver.h"

int main(int agrc, char** argv) {
	DiophantineSolver solver;


	solver.input(filestr(string(argv[1])));
	solver.solve();
	solver.output(filestr(string(argv[2])));
}