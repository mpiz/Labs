#include "DiophantineSolver.h"

int main(int agrc, char** argv) {
	DiophantineSolver solver;

	solver.input(string(argv[1]));
	solver.solve();
	solver.output(string(argv[2]));
}