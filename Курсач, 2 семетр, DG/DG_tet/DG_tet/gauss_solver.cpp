#include "gauss_solver.h"
#include <cstring>
#include <iostream>
#include <fstream>

using namespace std;

Gauss::Gauss()
{
    A = NULL;
}

Gauss::~Gauss()
{
    if(A)
    {
        for(int i = 0; i < n; i++)
            delete [] A[i];
        delete [] A;
    }
}

// Инициализация несимметричная
void Gauss::init(int * s_ig, int * s_jg,
                 double * s_gu, double * s_gl,
                 double * s_di, double * s_rp,
                 int s_n)
{
    n = s_n;

    A = new double * [n];
    for(int i = 0; i < n; i++)
    {
        A[i] = new double [n + 1];
        memset(A[i], 0, sizeof(double) * (n + 1));
    }

    for(int i = 0; i < n; i++)
    {
        for(int k = s_ig[i]; k < s_ig[i + 1]; k++)
        {
            A[i][s_jg[k]] = s_gl[k];
            A[s_jg[k]][i] = s_gu[k];
        }
        A[i][i] = s_di[i];
        A[i][n] = s_rp[i];
    }
}

// Инициализация симметричная
void Gauss::init(int * s_ig, int * s_jg,
                 double * s_di, double * s_gg,
                 double * s_rp, int s_n)
{
    n = s_n;

    A = new double * [n];
    for(int i = 0; i < n; i++)
    {
        A[i] = new double [n + 1];
        memset(A[i], 0, sizeof(double) * (n + 1));
    }

    for(int i = 0; i < n; i++)
    {
        for(int k = s_ig[i]; k < s_ig[i + 1]; k++)
        {
            A[i][s_jg[k]] = s_gg[k];
            A[s_jg[k]][i] = s_gg[k];
        }
        A[i][i] = s_di[i];
        A[i][n] = s_rp[i];
    }
}

// Получение решения
void Gauss::solve(double * solution)
{
    long N = n;
    //верхний треугольный вид
    for(long i = 0; i < N; i++)
    {
        if(!A[i][i])
        {
            bool flag = false;
            for(long j = i + 1; j < N && !flag; j++)
                if(A[j][i])
                {
                    for(long k = i; k <= N; k++)
                    {
                        double tmp = A[i][k];
                        A[i][k] = A[j][k];
                        A[j][k] = tmp;
                    }
                    flag = true;
                }
        }
        for(long j = N; j >= i; j--)
            A[i][j] = A[i][j] / A[i][i];
        for(long j = i + 1; j < N; j++)
            for(long k = N; k >= i; k--)
                A[j][k] -= A[i][k] * A[j][i];
    }
    //диагональный вид
    for(long i = N - 1; i > 0; i--)
        for(long j = i - 1; j >= 0; j--)
            A[j][N] -= A[j][i] * A[i][N];

    for(int i = 0; i < n; i++)
        solution[i] = A[i][n];
}
