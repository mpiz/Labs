#include "VFEM_E.h"

//#if defined _WIN32
//#include <windows.h>
//long mtime()
//{
//    return GetTickCount();
//}
//#else
//#include <ctime>
//long mtime()
//{
//    struct timespec t;
//    clock_gettime(CLOCK_MONOTONIC, & t);
//    return (long)t.tv_sec * 1000 + t.tv_nsec / 1000000;
//}
//#endif

const int triangle_el = 2;
const int tethedron_el = 4;

// == Основные функции ==

void VFEM_E::input_mesh(string inp_file) {
	ifstream inp_f(inp_file.c_str());

	if (!inp_f.is_open())
		throw;

	cout << "Skiping to nodes\n";
	string read_str;
	while(read_str != "$Nodes")
		inp_f >> read_str;

	//Считывание узлов
	inp_f >> nodes_n;

	cout << "Nodes memory: " << sizeof(node) * nodes_n << endl;
	node read_n;
	nodes.resize(nodes_n);

	cout << "Reading nodes. Total nodes: " << nodes_n << "\n";

	for(int i = 0; i < nodes_n; i++) {
		inp_f >> read_n;
		read_n.number--;
		nodes[i]  = read_n;
	}

	cout << "Skiping to elements\n";
	//Пролистываем
	while(read_str != "$Elements")
		inp_f >> read_str;

	//Считывание граней и элементов

	int total_els, tmp_int, el_type;
	inp_f >> total_els;
	int ph_area;
	int params_n;
	array<int, 3> trnode;
	array<int, 4> tetnode;

	cout << "Reading elements: " << total_els << endl << "Size of element: " << sizeof(tetelement) << endl;

	elements.reserve(total_els);
	
	for(int iter = 0; iter < total_els; iter++) {
		inp_f >> tmp_int >> el_type;
		inp_f >> params_n;
		inp_f >> ph_area;

		if((iter * 100 / total_els) % 5 == 0) {
			cout << iter << "\\" << total_els << "\r";
			cout.flush();
		}


		switch(el_type) {
		case triangle_el: {
			for(int i = 0; i < params_n-1; i++);
				inp_f >> tmp_int;
	
			for(int i = 0; i < 3; i++) {
				inp_f >> trnode[i];
				trnode[i]--;
			}
		} break;

		case tethedron_el: {
			for(int i = 0; i < params_n-1; i++);
				inp_f >> tmp_int;
	
			for(int i = 0; i < 4; i++) {
				inp_f >> tetnode[i];
				tetnode[i]--;
			}

			sort(tetnode.begin(), tetnode.end());

			elements.push_back(tetelement(nodes[tetnode[0]], nodes[tetnode[1]], nodes[tetnode[2]], nodes[tetnode[3]]));
			int add_num = elements.size() - 1;
		} break;

		}; 

		if (iter%10000 == 0 || iter == total_els-1) {
			cout << "\t" << iter << "\t\\" << total_els << "\n";
		}

	}
	cout << endl;

	elements_n = elements.size();

	//Построение дерева поиска

	ifstream gabs("gabarits.txt");
	double gx0, gx1, gy0, gy1, gz0, gz1;

	gabs >> gx0 >> gx1 >> gy0 >> gy1 >> gz0 >> gz1;

	gabs.close();

	cout << "Search tree construct begin" << endl;

	search_tree.construct_tree(gx0, gx1, gy0, gy1, gz0, gz1, elements);
	cout << "Search tree construct complete" << endl;
}

#pragma optimize( "", off )
tetelement* VFEM_E::function_in_point_tree(double x, double y, double z) {

/*	LARGE_INTEGER start, stop, timetime, fr;
	double time;
	QueryPerformanceFrequency(&fr);
	QueryPerformanceCounter(&start);*/

	tetelement* point_el = search_tree.find_point(x, y, z);

/*	QueryPerformanceCounter(&stop);
	timetime.QuadPart = stop.QuadPart - start.QuadPart;
	time = (double)timetime.QuadPart / (double)fr.QuadPart;*/
	return point_el;
}
#pragma optimize( "", on )

double VFEM_E::function_in_point_linear(double x, double y, double z) {

	LARGE_INTEGER start, stop, timetime, fr;
	double time;
	QueryPerformanceFrequency(&fr);
	QueryPerformanceCounter(&start);

	for (int el_i = 0; el_i < elements_n; el_i++) {
		if (elements[el_i].in_element(x,y,z))
			break;
	}
	QueryPerformanceCounter(&stop);
	timetime.QuadPart = stop.QuadPart - start.QuadPart;
	time = (double)timetime.QuadPart / (double)fr.QuadPart;

	return time;

}
