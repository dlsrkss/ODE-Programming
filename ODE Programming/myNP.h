/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : ¾ç µÎ ¿ø
StudentID        : 21900427
Created          : 26-03-2018
Modified         : 11-23-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.h
----------------------------------------------------------------*/

#ifndef		_MY_NP_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NP_H

#include "myMatrix.h"

extern double myfunc(const double t, const double y);
extern void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
extern void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0);
void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0, int method);
void sys2RK2(void func(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init);
extern void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P);
extern void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x);
extern void backsub(Matrix U, Matrix y, Matrix x);
extern void fwdsub(Matrix L, Matrix b, Matrix y);
Matrix newtonRoots(Matrix X0, Matrix myJacob(Matrix X), Matrix myFunc(Matrix X), float tol);
#endif