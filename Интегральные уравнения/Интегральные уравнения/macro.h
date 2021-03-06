#pragma once

#include <array>
#include <complex>
#include <vector>
#include <omp.h>
#include <functional>

using namespace std;
// ��� �������� ���������� ������ ���������� ����������
#define DEBUGOUTP

//������ ����� ���������� ������ ������� �������
//#define O2M2

// ��������� ��� ��������� ���������
#define GEOCONST 1e-6

//typedef complex<double> dcomplex;

typedef complex<double> dcomplex;

#ifndef matrix
#define matrix(n) array<array<double,(n)>,(n)>  //���������� ������� ������� n
#endif

#ifndef cmatrix
#define cmatrix(n) array<array<dcomplex,(n)>,(n)>  //���������� ������� ������� n
#endif

#ifndef recmatrix
#define recmatrix(n,m) array<array<double,(m)>,(n)>  //���������� ������� ������� n
#endif


typedef size_t dof_type;
typedef dof_type dof_info;
typedef vector<vector<dcomplex>> dyn_matrix;  // ������������ �������

#define geom_err 1e-6
