#ifndef _MATRIX_
#define _MATRIX_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace std;

class Matrix{
public:
	int n_row,n_column;
	double **data;
	
	//Constructores
	Matrix(const int v_size);
	Matrix(const int n_row,const int n_column);
	
	
	
	// Member Operators
	double& operator () (const int n);
	double& operator () (const int n_row,const int n_column);
	Matrix& operator + (Matrix &m);
	Matrix& operator + (double s);
	Matrix& operator - (Matrix &m);
	Matrix& operator - (double s);
	Matrix& operator * (Matrix &m);
	Matrix& operator * (double s);
	Matrix& operator / (Matrix &m);
	Matrix& operator / (double s);
	Matrix& operator = (Matrix &m);
	
	Matrix& inv();
	double det() const;

	
	
	// Non-member operators
	friend ostream& operator << (ostream &o,Matrix&m);
	
	
	
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);

#endif