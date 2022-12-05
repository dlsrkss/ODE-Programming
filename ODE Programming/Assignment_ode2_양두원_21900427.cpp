/*------------------------------------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University
Author          : 양 두 원
Student ID      : 21900427
Created         : 01-04-2019
Modified        : 11-23-2022
Language/ver	: C in MSVS2022
Course			: Numerical Programming 2022
Description     : Assignment_eigenvalue_ID.cpp
/------------------------------------------------------------------------------------------*/

#include "myNP.h"
#define PI 3.1415926
#define EU 0
#define RK2 1
#define RK3 2
Matrix myFunc(Matrix X);
Matrix myJacob(Matrix X);
Matrix myFunc2(Matrix X);
Matrix myJacob2(Matrix X);



Matrix myFunc(Matrix X)
{
	Matrix A = copyMat(X);
	A.at[0][0] = X.at[1][0] - 0.5 * (exp(X.at[0][0] / 2) + exp(-X.at[0][0] / 2));
	A.at[1][0] = 9 * X.at[0][0] * X.at[0][0] + 25 * X.at[1][0] * X.at[1][0] - 225;
	return A;
}

Matrix myJacob(Matrix X)
{

	Matrix A = zeros(X.rows, X.rows);
	A.at[0][0] = -0.25 * (exp(X.at[0][0] / 2) - exp(-X.at[0][0] / 2));
	A.at[0][1] = 1;
	A.at[1][0] = 18 * X.at[0][0];
	A.at[1][1] = 50 * X.at[1][0];
	return A;
}

Matrix myFunc2(Matrix X)
{
	
	Matrix A = copyMat(X); //출력 theta, delta x, delta y

		A.at[0][0] = cos(X.at[0][0]) * 0 - sin(X.at[0][0]) * 100 + X.at[1][0]-50;
		A.at[1][0] = sin(X.at[0][0]) * 0 + cos(X.at[0][0]) * 100 + X.at[2][0]-186.6025;
		A.at[2][0] = cos(X.at[0][0]) * 0 - sin(X.at[0][0]) * (-100) + X.at[1][0] - 150;
	
	return A;
}

Matrix myJacob2(Matrix X)
{
	Matrix A = createMat(X.rows, X.rows); //출력 theta, delta x, delta y

	A.at[0][0] = -sin(X.at[0][0]) * 0 - cos(X.at[0][0]) * 100;
	A.at[0][1] = 1;
	A.at[0][2] = 0;
	A.at[1][0] = cos(X.at[0][0]) * 0 - sin(X.at[0][0]) * 100;
	A.at[1][1] = 0;
	A.at[1][2] = 1;
	A.at[2][0] = -sin(X.at[0][0]) * 0 - cos(X.at[0][0]) * (-100);
	A.at[2][1] = 1;
	A.at[2][2] = 0;
	return A;
}

int main(int argc, char* argv[])
{
	/*==========================================================================*/
	/*					Variables declaration & initialization					*/
	/*--------------------------------------------------------------------------*/
	/*   - You can change the variable names									*/
	/*   - However, you must use the specified file name						*/
	/*	   : For each assignment, the file name will be notified on HISNET		*/
	/*==========================================================================*/
	double t0 = 0; double tf = 1;
	double h = 0.01;
	double N = (tf - t0) / h + 1;
	double t = 0;
	//double y1;
	//double y2[2] = {0 , 0};
	double* y1 = (double*)malloc(sizeof(double) * N); //y함수
	double* y2 = (double*)malloc(sizeof(double) * N); //z함수
	double y1_init = 0;
	double y2_init = 0.2;
	

	printf("-------------------------------------------------------------------------\n");
	printf("                               ODE	Results					             \n");
	printf("-------------------------------------------------------------------------\n");
	Matrix X1= zeros(2, 1);
	X1.at[0][0] = 2.5;
	X1.at[1][0] = 2;
	float tol = 1e-9;
	newtonRoots(X1, myJacob, myFunc, tol);

	Matrix X2 = zeros(3, 1);
	X2.at[0][0] = 0;
	X2.at[1][0] = 0;
	X2.at[2][0] = 0;
	newtonRoots(X2, myJacob2, myFunc2, tol);
	//sys2RK2(mckfunc, y1, y2, t0, tf, h, y1_init, y2_init);



}


void mckfunc(const double t, const double Y[], double dYdt[])
{
	double m = 1; double c = 7; double k = 6.9; double f = 5;
	double Fin = 2 * cos(2 * PI * f * t); 
	dYdt[0] = Y[1]; //1번째 풀어야 할 미방식


	// EXERCISE: MODIFY HERE
	dYdt[1] = -Y[1] * c / m - Y[0] * k / m + Fin / m; //z값
}


