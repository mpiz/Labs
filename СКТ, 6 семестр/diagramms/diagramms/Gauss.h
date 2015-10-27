#pragma once

#include <math.h>

bool SLAE_solution_Gauss(double **A, double *b, double *x, int n); //решение СЛАУ, метод Гаусса-Жордана

bool trianglematrix1(double **A, double *x, int n); //Приведение матрицы к верхнему треугольному виду
void transf1(double **A, double *x, int i,int n); // перестановка строк