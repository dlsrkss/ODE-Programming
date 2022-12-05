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
### Trapezode
_Integration using Trapezoidal method for discrete data inputs_
<img src="![캡처](https://user-images.githubusercontent.com/116098075/205694781-378e5915-388a-4c7c-8e41-ca265aacbded.PNG)
" width="300" height="200"> 
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
