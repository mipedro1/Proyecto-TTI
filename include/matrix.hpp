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
	Matrix();
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
	void transpose();
	
	double dot(Matrix& m);
	Matrix& cross(Matrix& m);
	Matrix& extract_vector(int indiceInicio, int indiceFinal);
	Matrix& union_vector(Matrix& m);
	Matrix& extract_row(int n);
	Matrix& extract_column(int n);
	Matrix& assign_row(int n, Matrix& m);
	Matrix& assign_column(int n, Matrix& m);
	
	
	// Non-member operators
	friend ostream& operator << (ostream &o,Matrix&m);
	
	
	
};

// Operator overloading
ostream& operator << (ostream &o, Matrix &m);

// Methods
Matrix& zeros(const int n_row, const int n_column);

Matrix& eye(int n);

Matrix& zeros(int n);

double norm(Matrix& m);


#endif