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

void mckfunc(const double t, const double Y[], double dYdt[]);


double myLossEx1(double x);
double myLossGradEx1(double x);

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);
double myLossEx2(Matrix Z);
Matrix myLossGradEx2(Matrix Z);



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
	
	/////////////
	int itrMax = 20;
	double	eta = 0.25;
	int k = 0;
	double grad_e = 0;
	double x = 0;
	double loss = 0;

	// Initial Condition
	x = 0;

	/************      Gradient Descent							     ************/

	for (k = 0; k < itrMax; k++) {
		grad_e = myLossGradEx1(x);

		// [TO-DO]  YOUR CODE GOES HERE
		// update x as  x=x-eta*gradL 
		x += -eta * grad_e;

		loss = myLossEx1(x);
		printf("iter =%d \t x=%0.3f \t L(x)=%0.3f\n", k, x, loss);
	}

	printf("\n\n");
	//////////////
	/************      Variables declaration & initialization      ************/

	itrMax = 20;
	eta = 0.25;
	Matrix grad_L;
	Matrix Z;
	Matrix H = zeros(2, 1);
	loss = 0;

	// Initial Condition
	double z0[] = { -2,-2 };
	Z = arr2Mat(z0, 2, 1);


	/************      Gradient Descent							     ************/

	for (k = 0; k < itrMax; k++) {
		grad_L = myLossGradEx2(Z);

		// [TO-DO]  YOUR CODE GOES HERE
		// update Z as  Z=Z+H, where H=-eta*gradL 
		
		H = mulScalar(grad_L, -eta);
		for (int i = 0; i < Z.rows; i++)
		{
			Z.at[i][0] = Z.at[i][0] + H.at[i][0];
		}

		loss = myLossEx2(Z);
		printf("iter =%d \t x=%0.3f \t y=%0.3f \t L(x)=%0.3f\n", k, Z.at[0][0], Z.at[1][0], loss);
	}

	printf("\n\n");
	/////////////////////////////////
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
	//newtonRoots(X2, myJacob2, myFunc2, tol);
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


double myLossEx1(double x) {
	double L = x * x - 4 * x + 6;
	return L;
}

double myLossGradEx1(double x) {

	double dL = 0;
	// [TO-DO] YOUR CODE GOES HERE
	dL = 2 * x - 4;


	return dL;
}

double myLossEx2(Matrix Z) {

	double L = 0;
	int n = Z.rows;
	double z1 = Z.at[0][0];
	double z2 = Z.at[1][0];

	// [TO-DO] YOUR CODE GOES HERE
	L = 3 * (z1 - 2) * (z1 - 2) + (z2 - 2) * (z2 - 2);


	return L;
}


Matrix myLossGradEx2(Matrix Z) {

	int n = Z.rows;
	Matrix dL = zeros(n, 1);

	double z1 = Z.at[0][0];
	double z2 = Z.at[1][0];

	// [TO-DO] YOUR CODE GOES HERE
	dL.at[0][0] = 6 * (z1 - 2);
	dL.at[1][0] = 2 * (z2 - 2);

	return dL;
};


Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}