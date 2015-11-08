#pragma once
#define _CRT_SECURE_NO_WARNINGS

#include "elements_classes.h"
#include "octal_tree.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

struct media {
	double mu;
	double epsilon;
	double sigma;

};

class VFEM_E {
 public:
	 void input_mesh(string inp_file);
	 double function_in_point_tree(double x, double y, double z);
	 double function_in_point_linear(double x, double y, double z);

 private:

	 int nodes_n;
	 int elements_n;


	 vector<node> nodes;
	 vector<tetelement> elements;

	 octal_tree<tetelement> search_tree;

};

