/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : 양 두 원
StudentID        : 21900427
Created          : 26-03-2018
Modified         : 11-23-2022
Language/ver     : C++ in MSVS2019

Description      : myNP.cpp
----------------------------------------------------------------*/

#include "myNP.h"

//EU
void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{

	int N = (tf - t0) / h + 1;
	double t = t0;
	double* yE;
	yE = (double*)malloc(sizeof(double) * N);
	yE[0] = y0;

	for (int i = 0; i < N; i++)
	{
		
		y[i] = yE[i];
		yE[i + 1] = yE[i] + myfunc(t, yE[i]) * h;
		printf("t[%d] : %f    yE[%d] : %f\n\n", i, t, i, yE[i]);
		t = t + h;
	}


}

void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{
	int N = (tf - t0) / h + 1;
	double t1 = t0;
	double t2 = 0;
	double* yEM;
	yEM = (double*)malloc(sizeof(double) * N);
	yEM[0] = y0;

	for (int i = 0; i < N; i++)
	{
		t2 = t1 + h;
		double slope1 = myfunc(t1, yEM[i]);
		double slope2 = myfunc(t2, yEM[i] + myfunc(t1, yEM[i])*h);
		yEM[i + 1] = yEM[i] + (slope1 + slope2) * h / 2;
		printf("t[%d] : %f    yEM[%d] : %f\n\n", i, t1, i, yEM[i]);
		t1 = t2;
	}

}
//RK2

void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{
	int N = (tf - t0) / h + 1;
	double t1 = t0;
	double t2 = 0;
	double alpha = 1;
	double beta = 1;
	double C1 = 0.5;
	double C2 = 0.5;
	double K1;
	double K2;
	
	double* yRK2;
	yRK2 = (double*)malloc(sizeof(double) * N);
	yRK2[0] = y0;

	for (int i = 0; i < N; i++)
	{
		t2 = t1 + h;
		K1 = myfunc(t1, yRK2[i]);
		K2 = myfunc(t1 + alpha * h, yRK2[i] + beta * K1 * h);
		yRK2[i + 1] = yRK2[i] + (C1 * K1 + C2 * K2) * h;
		printf("t[%d] : %f    yRK2[%d] : %f\n\n", i, t1, i, yRK2[i]);
		t1 = t2;
	}

}//RK2

void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
{
	int N = (tf - t0) / h + 1;
	double t1 = t0;
	double t2 = 0;
	double alpha2 = 0.5;
	double alpha3 = 1;
	double beta21 = 0.5;
	double beta31 = -1;
	double beta32 = 2;
	double C1 = 1 / 6;
	double C2 = 4 / 6;
	double C3 = 1 / 6;
	double K1;
	double K2;
	double K3;

	double* yRK3;
	yRK3 = (double*)malloc(sizeof(double) * N);
	yRK3[0] = y0;

	for (int i = 0; i < N; i++)
	{
		t2 = t1 + h;
		K1 = myfunc(t1, yRK3[i]);
		K2 = myfunc(t1 + alpha2 * h, yRK3[i] + beta21 * K1 * h);
		K3 = myfunc(t1 + alpha3 * h, yRK3[i] + beta31 * K1 * h + beta32 * K2 * h);
		yRK3[i + 1] = yRK3[i] + (K1 + 4 * K2 + K3) * h / 6;
		printf("t[%d] : %f    yRK3[%d] : %f\n\n", i, t1, i, yRK3[i]);
		t1 = t2;
	}

}//RK3

void ode(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0, int method)
{
	if (method == 0)
		odeEU(myfunc, y, t0, tf, h, y0);
	else if (method == 1)
		odeRK2(myfunc, y, t0, tf, h, y0);
	else if (method == 2)
		odeRK3(myfunc, y, t0, tf, h, y0);
}

void sys2RK2(void mckfunc(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
{
	//시간변수 정의
	int N = (tf - t0) / h + 1;
	double t1 = t0;
	double t2 = 0;

	//계수 정의
	double alpha = 1;
	double beta = 1;
	double gamma = 1;
	double C1 = 0.5;
	double C2 = 0.5;

	//행렬 선언
	double y[2];
	double z[2];
	double y_weight[2];
	double z_weight[2];

	//K선언
	double Ky1;
	double Kz1;
	double Ky2;
	double Kz2;

	y1[0] = y1_init;
	y2[0] = y2_init;

	printf("t[0] : 0    ysRK2[0] : %f\n\n", y1[0]);


	for (int i = 0; i < N - 1; i++)
	{
		t2 = t1 + h;
		y[0] = y1[i];
		y[1] = y2[i];
		
		mckfunc(t1, y, z);
		Ky1 = z[0];
		Kz1 = z[1];

		y_weight[0] = y1[i] + beta * Ky1 * h; // y의 가중치
		y_weight[1] = y2[i] + gamma * Kz1 * h; // y'의 가중치

		mckfunc(t1+alpha*h, y_weight, z_weight);
		Ky2 = z_weight[0];
		Kz2 = z_weight[1];

		y1[i + 1] = y1[i] + (C1 * Ky1 + C2 * Ky2) * h;
		y2[i + 1] = y2[i] + (C1 * Kz1 + C2 * Kz2) * h;
		t1 = t2;
		printf("t[%d] : %f    ysRK2[%d] : %f\n\n", i + 1, t1, i + 1, y1[i + 1]);
		
	}

}


//Matrix X0 = F(x) n*1
Matrix newtonRoots(Matrix X0, Matrix myJacob(Matrix X), Matrix myFunc(Matrix X), float tol)
{
	
	
	float error = 1e+5;
	double sum = 0;
	int Nmax = 100;
	int itcount = 0;
	double dx = 1e-9;
	Matrix X = copyMat(X0);
	Matrix J = myJacob(X0);
	Matrix F = myFunc(X0);
	
	Matrix H = zeros(J.rows, 1);
	Matrix U = zeros(J.rows, J.cols);
	Matrix L = zeros(J.rows, J.cols);
	Matrix P = eye(J.rows, J.cols);



	
	do {
		sum = 0;
	

		for (int i = 0; i < F.rows; i++)
		{
			F.at[i][0] = -F.at[i][0];
		}

		LUdecomp(J, L, U, P);
		solveLU(L, U, P, F, H);

		for(int i = 0; i <F.rows; i++)
		{
			X.at[i][0] = X.at[i][0] + H.at[i][0];
		}

		
		copyVal(myJacob(X), J);
		copyVal(myFunc(X), F);


		for (int i = 0; i < F.rows; i++)
		{
			sum += F.at[i][0] * F.at[i][0];
		}
		error = sum/F.rows;


		itcount++;
		
		printf("iteration:%d\t", itcount);
		printMat(X, "X");
		printf("Tolerance : %.10f\n", error);
	} while (itcount < Nmax && error > tol);

	
	return X;
}


void backsub(Matrix U, Matrix y, Matrix x)
{
	double s;
	for (int i = U.rows - 1; 0 <= i; i--) {
		s = 0;
		for (int j = i + 1; j < U.cols; j++)
		{
			s = s + U.at[i][j] * x.at[j][0];
		}
		x.at[i][0] = (y.at[i][0] - s) / (U.at[i][i]);
	}
}



void fwdsub(Matrix L, Matrix b, Matrix y)
{
	copyVal(b, y);
	double s = 0;
	for (int i = 0; i < L.cols; i++)
	{
		for (int j = i + 1; j < L.rows; j++)
		{
			s = L.at[j][i] * y.at[i][0];
			y.at[j][0] -= s;
		}
	}
}


void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P)
{
	double M;
	copyVal(zeros(L.rows, L.cols), L);
	copyVal(A, U);
	Matrix Temp = zeros(1, A.cols);
	int piv, max1, max2;
	double sp1 = 0; double sp2 = 0;
	//예외처리
	if (A.rows != A.cols) printf("ERROR : A is not N*N matrix");

	for (int k = 0; k < U.cols - 1; k++)
	{
		piv = k;
		max1 = k + 1;
		max2 = k + 1;

		//scaled pivoting
		for (int i = k; i < U.rows - 1; i++)
		{


			for (int j = k; j < U.cols; j++)
			{
				if (fabs(U.at[piv][max1]) < fabs(U.at[piv][j]))
					max1 = j;
			}
			sp1 = fabs(U.at[piv][k]) / fabs(U.at[piv][max1]);



			for (int j = k; j < U.cols; j++)
			{
				if (fabs(U.at[i + 1][max2]) < fabs(U.at[i + 1][j]))
					max2 = j;
			}
			sp2 = fabs(U.at[i + 1][k]) / fabs(U.at[i + 1][max2]);



			if (sp1 < sp2)
			{
				piv = i + 1;
				max1 = k + 1;
				max2 = k + 1;
			}
		}

		for (int j = 0; j < U.cols; j++)
		{
			double temp;
			temp = U.at[piv][j];
			U.at[piv][j] = U.at[k][j];
			U.at[k][j] = temp;

			double Ltemp;
			Ltemp = L.at[piv][j];
			L.at[piv][j] = L.at[k][j];
			L.at[k][j] = Ltemp;

			double ptemp;
			ptemp = P.at[piv][j];
			P.at[piv][j] = P.at[k][j];
			P.at[k][j] = ptemp;
		}

		//Gauss Elmination

		for (int i = k + 1; i < U.rows; i++)
		{
			//check for div by zero
			if (U.at[k][k] == 0) printf("ERROR : div value is zero");


			M = (U.at[i][k] / U.at[k][k]);
			L.at[i][k] = M;
			for (int j = k; j < U.cols; j++)
			{
				U.at[i][j] += U.at[k][j] * (-M);
			}
		}

	}
	//L에 identity Matirx add
	for (int i = 0; i < L.rows; i++)
	{
		L.at[i][i] += 1;
	}
}

void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x)
{

	Matrix y = zeros(L.rows, 1);
	fwdsub(L, b, y);
	backsub(U, y, x);
}