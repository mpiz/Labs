#pragma once

#include <array>
#include <complex>
#include <vector>
#include <omp.h>
#include <functional>

using namespace std;
// ƒл€ быстрого отключени€ вывода отладочной информации
#define DEBUGOUTP

//¬торой метод построение базиса второго пор€дка
//#define O2M2

//  онстанта дл€ сранвени€ геометрии
#define GEOCONST 1e-6

//typedef complex<double> dcomplex;

typedef complex<double> dcomplex;

#ifndef matrix
#define matrix(n) array<array<double,(n)>,(n)>  //квадратна€ матрица пор€дка n
#endif

#ifndef cmatrix
#define cmatrix(n) array<array<dcomplex,(n)>,(n)>  //квадратна€ матрица пор€дка n
#endif

#ifndef recmatrix
#define recmatrix(n,m) array<array<double,(m)>,(n)>  //квадратна€ матрица пор€дка n
#endif


typedef size_t dof_type;
typedef dof_type dof_info;
typedef vector<vector<dcomplex>> dyn_matrix;  // динамическа€ матрица

#define geom_err 1e-6
