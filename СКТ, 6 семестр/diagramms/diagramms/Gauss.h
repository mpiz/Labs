#pragma once

#include <math.h>

bool SLAE_solution_Gauss(double **A, double *b, double *x, int n); //������� ����, ����� ������-�������

bool trianglematrix1(double **A, double *x, int n); //���������� ������� � �������� ������������ ����
void transf1(double **A, double *x, int i,int n); // ������������ �����