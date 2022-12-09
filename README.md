# myMatrix Source code
>## createMat

```C
Matrix createMat(int _rows, int _cols);
```
>>- __parameter__
>>1. int __rows__ : rows of Matrix
>>2. int __cols__ : rows of Matrix

>__Example code__
```C
Matrix A = createMat(2, 3);
Matrix A = createMat(5, 5);
```


>## printMat

```C
void printMat(Matrix _A, const char* _name)
```
>>- __parameter__
>>1. ___A__ : printed Matrix
>>2. ___name__ : title

>__Example code__
```C
printMat(matA, "matA");
printMat(matB, "matB");
```


>## mulScalar

```C
Matrix mulScalar(Matrix _A, double B)
```
>>- __parameter__
>>1. ___A__ : subject Matrix A
>>2. ___B__ : scalar number

>__Example code__
```C
Matrix X = mulScalar(matA, 2);
Matrix X = mulScalar(matB, 5);
```


>## mulMat

```C
Matrix mulMat(Matrix _A, Matrix _B);
```
>>- __parameter__
>>1. __A__ : first Matrix 
>>2. __B__ : second Matrix

>__Example code__
```C
mulMat(matA, matB);
```

>## initMat

```C
void	initMat(Matrix _A, double _val);
```
>>- __parameter__
>>1. __A__ : subject Matrix 
>>2. __B__ : initial value

>__Example code__
```C
initMat(matA, 1);
```


>## zeros

```C
Matrix	zeros(int _rows, int _cols);
```
>>- __parameter__
>>1. ___rows__ : rows of Matrix
>>2. ___cols__ : columns of Matrix

>__Example code__
```C
Matrix X = zeros(3, 4);
```

>## ones
create Matrix of all zeros

```C
Matrix	ones(int _rows, int _cols);
```
>>- __parameter__
>>1. ___rows__ : rows of Matrix
>>2. ___cols__ : columns of Matrix

>__Example code__
```C
Matrix X = ones(3, 3);
```


>## transpose

```C
Matrix transpose(Matrix _A);
```
>>- __parameter__
>>1. __A__ : subject Matrix 

>__Example code__
```C
transpose(matA);
```

>## copyMat

```C
Matrix	copyMat(Matrix _A);
```
>>- __parameter__
>>1. __A__ : subject Matrix 

>__Example code__
```C
matB = copyMat(matA);
```

>## copyVal

```C
void copyVal(Matrix _A, Matrix _B);
```
>>- __parameter__
>>1. ___A__ : subject Matrix 
>>2. ___B__ : replicated Matrix 

>__Example code__
```C
copyVal(matA, matB);
```




# Practice
These functions are useful tools for making any function in Numerical Programming Class

### Numerical Programming API
```C
#include "myNP.h"
```

## Non-Linear Solver
### bisection()
1. select a, b
2. if f(a)*f(b) appear negative sign, there is a value
3. it selects (a+b)/2, a | if f((a+b)/2)*f(a)< 0 : b =x; 
```C
double bisection(double fn(double _x), double _a, double _b, double _tau)
```

#### parameter
* __fn__ : function
* __a__ : initial value 1
* __b__ : initial value 2
* __tau__ : error value

#### example code
```C
double init_1 = -3;
double init_2 = 5;
double tau = 0.00001;

BS_result = bisection(function, init_1, init_2, tau);
```



- - -
### NewtonRaphson()
1. select a
2. it finds roots as x-intercept of the straight line tangent to the curve at a 
```C
double newtonRaphson(double fn(double x), double dfn(double x), double _x0, double _tol)
```

#### parameter
* __fn__ : function
* __dfn__ : differential function
* __x0__ : initial value
* __tol__ : tolerance error

#### example code
```C
float tol = 0.00001;
double NR_result;
double x0 = 1;

NR_result = newtonRaphson(fn, dfn, x0, tol);
```
- - -



### Secant()

```C
double secant(double function(double _x), double _x0, double _x1, double _tau)
```

#### parameter
* __function__ : function
* __x0__ : initial value 1
* __x1__ : initial value 2
* __tau__ : error value

#### example code
```C
float tau = 0.00001;
double SC_result;
double x0 = 1;
double x1 = 3;

SC_result = secant(fn, x0, x1, tau);
```
- - -


### NewtonRoots()

```C
Matrix newtonRoots(Matrix X0, Matrix myJacob(Matrix X), Matrix myFunc(Matrix X), float tol)
```

#### parameter
* __X0__ : initial value Matrix
* __myJacob__ : Jacobian Matrix
* __myFunc__ : Function Matrix
* __tol__ : tolerance error

#### example code
```C
Matrix X= zeros(2, 1);
X.at[0][0] = 2.5;
X.at[1][0] = 2;
float tol = 1e-9;

newtonRoots(X, myJacob, myFunc, tol);
```
- - -


## Numerical Differentiation
### gradient1D()
_Differential operation based on data array_

```C
void gradient1D(double x[], double y[], double dydx[], int m) 
```

#### parameter
* __x[]__ : x data
* __y[]__ : y data
* __dydx[]__ : differential coefficient array to return
* __m__ : array 

#### example code
```C
int m = 21;
double t[21] = { 0 };
for (int i = 0; i < m; i++) t[i] = 0.2 * i;
double x[] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59 };
double  dxdt[21] = { 0 };

gradient1D(t, x, dxdt, m); //m : 데이터 갯수 dxdt : 업데이트할 변수
printVec(dxdt, m);
```
- - -



### gradientFunc()
_Differential operation based on function_
```C
void gradientFunc(double func(const double x), double x[], double dydx[], int m)
```

#### parameter
* __func__ : function
* __x[]__ : x data
* __dydx[]__ : differential coefficient array to return
* __m__ : array size

#### example code
```C
double xin = -5.87;
double y = myFunc(xin);
double dydx[21];

gradientFunc(myFunc, t, dydx, m);
printVec(dydx, m);
```
- - -

### acceleration()
_double Derivative operation based on function_
```C
void acceleration(double x[], double y[], double d2ydx2[], int m)
```

#### parameter
* __x[]__ : x data
* __y[]__ : y data
* __dy2dx2[]__ : double derivative coefficient array to return
* __m__ : array size

#### example code
```C
int m = 21;
double t[21] = { 0 };
for (int i = 0; i < m; i++) t[i] = 0.2 * i;
double x[] = { -5.87, -4.23, -2.55, -0.89, 0.67, 2.09, 3.31, 4.31, 5.06, 5.55, 5.78, 5.77, 5.52, 5.08, 4.46, 3.72, 2.88, 2.00, 1.10, 0.23, -0.59 };
double  dx2dt2[21] = { 0 };

acceleration(t, x, dx2dt2, m);
printVec(dx2dt2, m);
```
- - -






## Integarion
### IntegrateRect()
_Integration using Rectangular method for discrete data inputs_  
<img src="https://user-images.githubusercontent.com/116098075/205697135-40bdff8d-8cc6-40e9-9e24-6bc339c95ee8.png" width ="400" height="300"> 
```C
double IntegrateRect(double x[], double y[], int m)
```

#### parameter
* __x[]__ : x data
* __y[]__ : y data
* __m__ : array size

#### example code
```C
double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
int M = sizeof(x) / sizeof(x[0]);

double I_rect = IntegrateRect(x, y, M);
printf("I_rect = %f\n\n", I_rect);
```
- - -



### Trapezode()
_Integration using Trapezoidal method for discrete data inputs_
<img src="https://user-images.githubusercontent.com/116098075/205694781-378e5915-388a-4c7c-8e41-ca265aacbded.PNG" width="600" height="300"> 
```C
double trapz(double x[], double y[], int m)
```

#### parameter
* __x[]__ : x data
* __y[]__ : y data
* __m__ : array size

#### example code
```C
double x[] = { 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60 };
double y[] = { 0, 3, 8, 20, 33, 42, 40, 48, 60, 12, 8, 4, 3 };
int M = sizeof(x) / sizeof(x[0]);

double I_trapz = trapz(x, y, M);
printf("I_trapz = %f\n\n", I_trapz);
```
- - -




### Simpson13()
_Integration using Simpson 1/3 method for discrete data inputs_  
<img src="https://user-images.githubusercontent.com/116098075/205697964-9efb339c-109c-4c02-8bea-212a5b27f5de.png" width="400" height="300"> 
```C
double simpson13(double x[], double y[], int m) {
```

#### parameter
* __x[]__ : x data
* __y[]__ : y data
* __m__ : array size

#### example code
```C
double x[] = { -3, -2.25, -1.5, -0.75, 0, 0.75, 1.5, 2.25, 3 };
double y[] = { 0, 2.1875, 3.75, 4.6875, 5, 4.6875, 3.75, 2.1875, 0 };
int M = sizeof(x) / sizeof(x[0]);

double I_simpson = simpson13(x, y, M);
	printf("I_simpson13  = %f\n\n", I_simpson);
```
- - -



### Simpson38()
_Integration using Simpson 3/8 method for discrete data inputs_  
<img src="https://user-images.githubusercontent.com/116098075/205698721-d2322603-efe3-404c-946a-8872bc25e084.png" width="400" height="300"> 
```C
double simpson38(double x[], double y[], int m) {
```

#### parameter
* __x[]__ : x data
* __y[]__ : y data
* __m__ : array size

#### example code
```C
double x[] = { -3, -2.25, -1.5, -0.75, 0, 0.75, 1.5, 2.25, 3 };
double y[] = { 0, 2.1875, 3.75, 4.6875, 5, 4.6875, 3.75, 2.1875, 0 };
int M = sizeof(x) / sizeof(x[0]);

double I_simpson = simpson38(x, y, M);
	printf("I_simpson38  = %f\n\n", I_simpson);
```
- - -



### Function integration()
_Integration using Simpson 1/3 method for function inputs_
```C
double integral(double func(const double _x), double a, double b, int n)
```

#### parameter
* __func__ : function
* __a__ : initial value 1
* __b__ : initial value 2
* __n__ : array size

#### example code
```C
int a = -1;
int b = 1;
int n = 12;

double Integrals = integral(myFunc, a, b, n);
printf("func_simpson13  = %f\n\n", Integrals);
```
- - -



## Linear Solver
### gaussElim()
_solves for vector x from Ax=b, a linear system problem Using Gauss Elimination_

```C
gaussElim(Matrix A, Matrix vecb, Matrix U, Matrix d, Matrix P)
```

#### parameter
* __A__ : Matrix A in structure Matrix form. Should be (nxn) square
* __vecb__ : vector vecb in structure Matrix form. Should be (nx1)
* __U__ :  Matrix U in structure Matrix form. Should be (nxn) square
* __d__ : vector d in structure Matrix form. Should be (nx1)
* __P__ : Partial Pivoting P in structrue Matrix form. Should be (nxn)

#### example code
_Matirx A and vecb data is already filled in notepad_
```C
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = copyMat(Q1_matA);
Matrix d = zeros(vecb.rows, 1);
Matrix P = eye(matA.rows, matA.cols);

gaussElim(A, vecb, U, d, P);
```
- - -



### LUdecomp()
_factors a matrix as the product of a lower triangular matrix and upper triangular matrix_

```C
void LUdecomp(Matrix A, Matrix L, Matrix U, Matrix P)
```

#### parameter
* __A__ : Matrix A in structure Matrix form. Should be (nxn) square
* __L__ :  Matrix U in structure Matrix form. Should be (nxn) square
* __U__ :  Matrix U in structure Matrix form. Should be (nxn) square
* __P__ : Partial Pivoting P in structrue Matrix form. Should be (nxn)

#### example code
_Matirx A and vecb data is already filled in notepad_
```C
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix matU = copyMat(Q1_matA);
Matrix matL = createMat(Q1_matA.rows, Q1_matA.cols);
Matrix P = eye(matA.rows, matA.cols);

LUdecomp(matA, matL, matU, P);
```
- - -



### fwdsub()
_find y using forward substitution_

```C
void fwdsub(Matrix L, Matrix b, Matrix y)
```

#### parameter
* __b__ :  vector b in structure Matrix form. Should be (nx1)
* __y__ :  vector y in structure Matrix form. Should be (nx1)
* __L__ :  Matrix L in structure Matrix form. Should be (nxn) square

#### example code
exsit LUdecomp before
```C
fwdsub(L, vecb, y);
```
- - -



### backsub()
_find x using back substitution_

```C
void backsub(Matrix U, Matrix y, Matrix x)
```

#### parameter
* __x__ :  vector x in structure Matrix form. Should be (nx1)
* __y__ :  vector y in structure Matrix form. Should be (nx1)
* __U__ :  Matrix U in structure Matrix form. Should be (nxn) square

#### example code
exsit LUdecomp before
```C
backsub(matU, y, x);
```
- - -



### solveLU()
_find x forward substitution, backward substitution_

```C
void solveLU(Matrix L, Matrix U, Matrix P, Matrix b, Matrix x)
```

#### parameter
* __L__ :  Matrix U in structure Matrix form. Should be (nxn) square
* __U__ :  Matrix U in structure Matrix form. Should be (nxn) square
* __P__ : Partial Pivoting P in structrue Matrix form. Should be (nxn)
* __b__ :  vector b in structure Matrix form. Should be (nx1)
* __x__ :  vector x in structure Matrix form. Should be (nx1)

#### example code
exsit LUdecomp before
```C
Matrix matL = createMat(matA.rows, matA.cols);
Matrix matU = copyMat(matA);
Matrix P = eye(Q1_matA.rows, Q1_matA.cols);
Matrix vecb = txt2Mat(path, "prob1_vecb");
Matrix x = zeros(vecb.rows, 1);

void solveLU(matL, matU, P, vecb, x);
```
- - -



### invMat()
_inverse Matrix_

```C
double invMat(Matrix A, Matrix Ainv)
```

#### parameter
* __A__ :  Matrix A in structure Matrix form. Should be (nxn) square
* __Ainv__ : Matrix Ainv in structure Matrix form. Should be (nxn) square

#### example code
invMat using LUdecomp, forward substitution, backward substitution
```C
Matrix matA = txt2Mat(path, "prob1_matA");
Matrix matAinv = copyMat(matA);

invMat(matA, matAinv);
```
- - -


## Eigenvalue & Eigenvector
### QRHousehold()
_QR factorization_

```C
void QRHousehold(Matrix A, Matrix Q, Matrix R)
```

#### parameter
* __A__ :  Matrix A in structure Matrix form. Should be (nxn) square
* __Q__ :  Matrix Q in structure Matrix form. Should be (nxn) square
* __R__ :  Matrix R in structure Matrix form. Should be (nxn) square

#### example code
```C
Matrix matA = txt2Mat(path, "prob_matA");
Matrix matQ = copyMat(matA);
Matrix matR = copyMat(matA);

QRHousehold(matA, matQ, matR);
```
- - -



### eig()
_eigenvalue_

```C
Matrix eig(Matrix A)
```

#### parameter
* __A__ :  Matrix A in structure Matrix form. Should be (nxn) square


#### example code
apply QRHousehold internal
```C
Matrix matA = txt2Mat(path, "prob_matA");

eigenvalue = eig(matA);
```
- - -




### eigvec()
_eigenvector_

```C
Matrix eigvec(Matrix A)
```

#### parameter
* __A__ :  Matrix A in structure Matrix form. Should be (nxn) square


#### example code
apply invMat internal
```C
Matrix matA = txt2Mat(path, "prob_matA");

eigenvector = eigvec(matA);
```
- - -


## ODE-IVP
### odeEU()
_solve ordinary differential equation using Euler explicit method ordinary_

```C
void odeEU(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
```
### odeEM()
_solve ordinary differential equation using Euler Modified method ordinary_

```C
void odeEM(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
```

### odeRK2()
_solve ordinary differential equation using Runge-Kutta second method ordinary_

```C
void odeRK2(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
```

### odeRK3()
_solve ordinary differential equation using Runge-Kutta third method ordinary_

```C
void odeRK3(double myfunc(const double t, const double y), double y[], double t0, double tf, double h, double y0)
```
#### parameter
* __myfunc__ :  function with t, y
* __y[]__ :  allocated array y
* __t0__ : initial time value
* __tf__ : final time value
* __h__ : interval value
* __y0__ : initial y value


#### example code
```C
double  t0 = 0; double tf = 0.1;
double h = 0.001;
double N = (tf - t0) / h + 1;
double t = 0;
double* y = (double*)malloc(sizeof(double) * N);
double y0 = 0;

/*double myfunc(const double t, const double y)
{
	int i = 0;
	double tau = 1;
	double T = 1 / tau;
	double f = 10;
	double Vm = 1;
	double w = 2 * PI * f * t;
	double funcx = -T * y + 1 * Vm * cos(w);
	return funcx;
}*/

odeEU(myfunc, y, t0, tf, h, y0);
odeEM(myfunc, y, t0, tf, h, y0);
odeRK2(myfunc, y, t0, tf, h, y0);
odeRK3(myfunc, y, t0, tf, h, y0);
```
- - -



### sys2RK2()
_solve secondary differential equation using Runge-Kutta second method ordinary_  
<img src="https://user-images.githubusercontent.com/116098075/205846940-038b878e-ff40-4808-8d0b-46d4c32e89ee.png" width="300" height="80"> 
<img src="https://user-images.githubusercontent.com/116098075/205846677-8e1aabfa-f54f-4ad6-a2b8-5907786599f1.png" width="400" height="80"> 
<img src="https://user-images.githubusercontent.com/116098075/205847104-bf7a0b49-4287-45b5-a76a-6aae576f397e.png" width="400" height="150"> 
```C
void sys2RK2(void mckfunc(const double t, const double Y[], double dYdt[]), double y1[], double y2[], double t0, double tf, double h, double y1_init, double y2_init)
```

#### parameter
* __mckfunc__ :  function with t, y, dydt
* __y1[]__ :  y, ex) y[0] = y, y[1] = dy/dt
* __y2[]__ :  z, ex) z[0] = dy/dt, z[1] = d2y/dt2
* __t0__ : initial time value
* __tf__ : final time value
* __h__ : interval value
* __y1__ : initial y value
* __y2__ : initial dy/dt value


#### example code
```C
double t0 = 0; double tf = 1;
double h = 0.01;
double N = (tf - t0) / h + 1;
double t = 0;
double* y1 = (double*)malloc(sizeof(double) * N); //y함수
double* y2 = (double*)malloc(sizeof(double) * N); //z함수
double y1_init = 0;
double y2_init = 0.2;

/*
void mckfunc(const double t, const double Y[], double dYdt[])
{
	double m = 1; double c = 7; double k = 6.9; double f = 5;
	double Fin = 2 * cos(2 * PI * f * t); 
	dYdt[0] = Y[1]; //1번째 풀어야 할 미방식


	// EXERCISE: MODIFY HERE
	dYdt[1] = -Y[1] * c / m - Y[0] * k / m + Fin / m; //z값
}
*/

sys2RK2(mckfunc, y1, y2, t0, tf, h, y1_init, y2_init);
```
- - -



## Linear Regression
### linearRegression()
_straight line that fits the data_
```C
Matrix	linearRegression(Matrix _x, Matrix _y)
```

#### parameter
* __x__ :  x data array Matrix. should be (nx1)
* __y__ :  y data array Matrix. should be (nx1)
#### example code
```C
int M = 6; //number of data
double T_array[] = { 30, 40, 50, 60, 70, 80 };
double P_array[] = { 1.05, 1.07, 1.09, 1.14, 1.17, 1.21 };
Matrix T = arr2Mat(T_array, M, 1);
Matrix P = arr2Mat(P_array, M, 1);

Matrix z = linearRegression(T, P);
double answer = z.at[0][0] * 100 + z.at[1][0]; //particular y-point at x-point

printMat(T, "T");
printMat(P, "P");
printMat(z, "z");
printf("answer = %lf \n", answer);
```
- - -



### polyfit()
_curve that fits the data_
```C
void polyfit(Matrix _x, Matrix _y, Matrix _z, int n)
```

#### parameter
* __x__ :  x data array Matrix. should be (nx1)
* __y__ :  y data array Matrix. should be (nx1)
* __z__ : coefficient Matrix. should be ((n+1)x1)
* __n__ : curve function order
#### example code
```C
int n = 4; // curve function order
double strain[] = { 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0 };
double stress[] = {0, 3, 4.5, 5.8, 5.9, 5.8, 6.2, 7.4, 9.6, 15.6, 20.7, 26.7, 31.1, 35.6, 39.3, 41.5};

Matrix Strain = arr2Mat(strain, N, 1);
Matrix Stress = arr2Mat(stress, N, 1);
Matrix a = zeros(n+1, 1); //function coefficient

polyfit(Strain, Stress, a, n);
printMat(a, "a");
```
- - -


