#ifndef GAUSS_H_INCLUDED
#define GAUSS_H_INCLUDED

#include <cmath>
#include <cstdlib>

class Gauss
{
public:
    void init(int * s_ig, int * s_jg, // Инициализация несимметричная
              double * s_gu, double * s_gl,
              double * s_di, double * s_rp,
              int s_n);
    void init(int * s_ig, int * s_jg, // Инициализация симметричная
              double * s_di, double * s_gg,
              double * s_rp, int s_n);
    void solve(double * solution);          // Получение решения
    Gauss();
    ~Gauss();
private:
    int n;
    double ** A;
};


#endif // GAUSS_H_INCLUDED
