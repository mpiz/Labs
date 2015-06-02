#pragma once

#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <cmath>

class BCG_LU
{
private:
    double *Mdi,*Mggl,*Mggu;
    double *Di,*Ggl,*Ggu,*F,*Res;
    int *Ig,*Jg;
    double *Res0;
    int N,N_ggu;
    long int maxiter,iter;
    double eps;
    double nev_r,nev_r_;

public:
    //инициализация
    void init(int * s_ig, int * s_jg, double * s_gu, double * s_gl, double * s_di, double * s_rp, int s_n);
    void solve(double * solution); //Получение решения

    int flag;
    BCG_LU();
    ~BCG_LU();
    //void input();
    void bsg_LU();
    void output();

    void LU();
    void slau_U(double *x, double *f);
    void slau_Ut(double *x, double *f);
    void slau_L(double *x, double *f);
    void slau_Lt(double *x, double *f);

    void MultVM(double *x,double *y);
    void MultVM_T(double *x,double *y);
    void copy(double *x,double *y,int n);
    double norma(double *v);
    double scalmult(double *x,double *y);
};

