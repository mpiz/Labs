#pragma once

#include "macro.h"
#include <math.h>

bool SLAE_solution_Gauss(dyn_matrix& A, dcomplex *b, dcomplex *x, int n); //������� ����, ����� ������-�������

bool trianglematrix1(dyn_matrix& A, dcomplex *x, int n); //���������� ������� � �������� ������������ ����
void transf1(dyn_matrix& A, dcomplex *x, int i,int n); // ������������ �����