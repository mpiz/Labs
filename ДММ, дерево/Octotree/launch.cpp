#include "VFEM_E.h"
#include <stdlib.h>
#include <stdio.h>


int main() {


	VFEM_E octtree;
	octtree.input_mesh("web_spl2");

	int samples_n = 1e6;
	//double x, y, z;
	double tree_time = 0, linear_time = 0;
	double tt_iter, lt_iter;

	ofstream outp_info("iters_info.txt");

	srand(GetTickCount());

	vector<double> x, y, z;
	x.resize(samples_n);
	y.resize(samples_n);
	z.resize(samples_n);
	for(int i = 0; i < samples_n; i++) {
		x[i] = rand()%100 * 1.0 / 100;
		y[i] = rand()%100 * 1.0 / 100;
		z[i] = rand()%100 * 1.0 / 100;
	}

	LARGE_INTEGER start, stop, timetime, fr;
	QueryPerformanceFrequency(&fr);
	QueryPerformanceCounter(&start);

//	getchar();


	for (int i = 0; i < samples_n; i++) {
		auto el = octtree.function_in_point_tree(x[i], y[i], z[i]);
	}


	QueryPerformanceCounter(&stop);
	timetime.QuadPart = stop.QuadPart - start.QuadPart;
	tree_time = (double)timetime.QuadPart / (double)fr.QuadPart;

	printf("Avg time, tree method: %.3e, Full time: %.3e\n", tree_time / samples_n, tree_time);
	//printf("Avg time, linear method: %.3e, Full time: %.3e\n", linear_time / samples_n, linear_time);
	outp_info.close();

	system("pause");

	return 0;
}