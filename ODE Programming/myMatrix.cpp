/*----------------------------------------------------------------\
@ Numerical Programming by Young-Keun Kim - Handong Global University

Author           : �� �� ��
StudentID        : 21900427
Created          : 26-03-2018
Modified         : 11-23-2022
Language/ver     : C++ in MSVS2019

Description      : myMatrix.cpp
----------------------------------------------------------------*/

#include "myMatrix.h"


// Create Matrix with specified size
Matrix	createMat(int _rows, int _cols)
{
	// check matrix dimension
	if (_rows < 0 || _cols < 0) {
		printf("\n****************************************************");
		printf("\n  ERROR!!: dimension error at 'createMat' function");
		printf("\n****************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out;
	// 1. Allocate row array first
	Out.at = (double**)malloc(sizeof(double*) * _rows);
	// 2. Then, allocate column 
	for (int i = 0; i < _rows; i++)
		Out.at[i] = (double*)malloc(sizeof(double) * _cols);
	// 3. Initialize row & column values of a matrix
	Out.rows = _rows;
	Out.cols = _cols;

	return Out;
}

// Free a memory allocated matrix
void	freeMat(Matrix _A)
{
	// 1. Free allocated column memory
	for (int i = 0; i < _A.rows; i++)
		free(_A.at[i]);
	// 2. Free allocated row memory
	free(_A.at);
}

// Create a matrix from a text file
Matrix	txt2Mat(std::string _filePath, std::string _fileName)
{
	std::ifstream file;
	std::string temp_string, objFile = _filePath + _fileName + ".txt";
	int temp_int = 0, nRows = 0;

	file.open(objFile);
	if (!file.is_open()) {
		printf("\n*********************************************");
		printf("\n  Could not access file: 'txt2Mat' function");
		printf("\n*********************************************\n");
		return createMat(0, 0);
	}
	while (getline(file, temp_string, '\t'))
		temp_int++;
	file.close();

	file.open(objFile);
	while (getline(file, temp_string, '\n'))
		nRows++;
	file.close();

	int nCols = (temp_int - 1) / nRows + 1;
	Matrix Out = createMat(nRows, nCols);

	file.open(objFile);
	for (int i = 0; i < nRows; i++)
		for (int j = 0; j < nCols; j++) {
			file >> temp_string;
			Out.at[i][j] = stof(temp_string);
		}
	file.close();

	return Out;
}

// Print matrix
void	printMat(Matrix _A, const char* _name)
{
	printf("%s =\n", _name);
	for (int i = 0; i < _A.rows; i++) {
		for (int j = 0; j < _A.cols; j++)
			printf("%15.6f\t", _A.at[i][j]);
		printf("\n");
	}
	printf("\n");
}

Matrix mulScalar(Matrix _A, double B)
{
	Matrix Out = zeros(_A.rows, _A.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
			Out.at[i][j] += B * _A.at[i][j];
	return Out;
}


Matrix mulMat(Matrix _A, Matrix _B)
{
	if (_A.cols != _B.rows)
	{
		printf("\n******************************************************");
		printf("\n  ERROR!!: dimension error at 'mulMat' function");
		printf("\n******************************************************\n");
		return createMat(0, 0);
	}

	Matrix Out = zeros(_A.rows, _B.cols);

	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _B.cols; j++)
			for (int k = 0; k < _B.rows; k++)
				Out.at[i][j] += _A.at[i][k] * _B.at[k][j];

	return Out;
}
// initialization of Matrix elements
void	initMat(Matrix _A, double _val)
{
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
		{
			_A.at[i][j] = _val;
		}
}

// Create matrix of all zeros
Matrix	zeros(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	initMat(Out, 0);

	return Out;
}

// Create matrix of all ones
Matrix	ones(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	initMat(Out, 1);

	return Out;
}

// Create identity 
Matrix	eye(int _rows, int _cols)
{
	Matrix Out = createMat(_rows, _cols);
	initMat(Out, 0);
	for (int i = 0; i < Out.rows; i++)
	{
		Out.at[i][i] = 1;
	}

	return Out;
}

// Create Transpose matrix
//extern	Matrix	transpose(Matrix _A);
Matrix transpose(Matrix _A)
{
	Matrix Out = createMat(_A.cols, _A.rows);
	for (int i = 0; i < _A.rows; i++)

		for (int j = 0; j < _A.cols; j++)
		{
			Out.at[j][i] = _A.at[i][j];

		}
	return Out;
}
// Copy matrix
Matrix	copyMat(Matrix _A)
{
	Matrix Out = createMat(_A.rows, _A.cols);
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
		{
			Out.at[i][j] = _A.at[i][j];
		}
	return Out;
}

Matrix augMat(Matrix _A, Matrix vectorb)
{
	Matrix Out = createMat(_A.rows, _A.cols + vectorb.cols);
	for (int i = 0; i < _A.rows; i++)
	{
		for (int j = 0; j < _A.cols + vectorb.cols; j++)
		{
			if (j < _A.cols)
				Out.at[i][j] = _A.at[i][j];
			else Out.at[i][j] = vectorb.at[i][0];
		}
	}
	return Out;
}

void copyVal(Matrix _A, Matrix _B)
{
	for (int i = 0; i < _A.rows; i++)
		for (int j = 0; j < _A.cols; j++)
		{
			_B.at[i][j] = _A.at[i][j];
		}
}

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}
// Copy matrix Elements from A to B
//extern	void	copyVal(Matrix _A, Matrix _B);