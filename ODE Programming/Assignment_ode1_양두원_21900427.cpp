/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University
Author          : ¾ç µÎ ¿ø
Student ID      : 21900427
Created         : 01-04-2019
Modified        : 11-15-2022
Language/ver	: C in MSVS2022
Course			: Numerical Programming 2022
Description     : Assignment_eigenvalue_ID.cpp
/------------------------------------------------------------------------------------------*/

double myfunc(const double t, const double y);
void mckfunc(const double t, const double Y[], double dYdt[]);
#include "myNP.h"
#define PI 3.1415926
#define EU 0
#define RK2 1
#define RK3 2

int main(int argc, char* argv[])
{
	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/
	double  t0 = 0; double tf = 0.1;
	double h = 0.001;
	double N = (tf - t0) / h + 1;
	double t = 0;
	double* y = (double*)malloc(sizeof(double) * N);
	double y0 = 0;
	


	printf("-------------------------------------------------------------------------\n");
	printf("                               ODE	Results					             \n");
	printf("-------------------------------------------------------------------------\n");

	odeEU(myfunc, y, t0, tf, h, y0);
	odeEM(myfunc, y, t0, tf, h, y0);
	odeRK2(myfunc, y, t0, tf, h, y0);
	odeRK3(myfunc, y, t0, tf, h, y0);
	//ode(myfunc, y, t0, tf, h, y0, EU);


}

double myfunc(const double t, const double y)
{
	int i = 0;
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double w = 2 * PI * f * t;
	double funcx = -T * y + 1 * Vm * cos(w);
	return funcx;
}


void mckfunc(const double t, const double Y[], double dYdt[])
{
	double m = 1; 
	double c = 7; 
	double k = 6.9; 
	double f = 5;
	double Fin = 2 * cos(2 * PI * f * t); 
	dYdt[0] = Y[1];

	// EXERCISE: MODIFY HERE
	dYdt[1] = ;
}