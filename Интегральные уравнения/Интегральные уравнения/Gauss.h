#pragma once

#include "macro.h"
#include <math.h>

bool SLAE_solution_Gauss(dyn_matrix& A, dcomplex *b, dcomplex *x, int n); //решение СЛАУ, метод Гаусса-Жордана

bool trianglematrix1(dyn_matrix& A, dcomplex *x, int n); //Приведение матрицы к верхнему треугольному виду
void transf1(dyn_matrix& A, dcomplex *x, int i,int n); // перестановка строк