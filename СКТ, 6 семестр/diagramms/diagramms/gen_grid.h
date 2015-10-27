#pragma once

#include <stdio.h>
#include <conio.h>
#include "math.h"
#include <vector>
#include <algorithm>
using namespace std;


void generate_1reg_grid(double a, double b, double h, double *grid_mass, int n);
void generate_1unreg_grid(double a, double b, double h, double kr, vector<double>&grid_mass, int n, bool flag);
bool generate_2reg_grid();
bool generate_2unreg_grid();

