#include "../include/matrix.hpp"


Matrix::Matrix(const int v_size){
	if(v_size<0){
		cout << "Vector create: error in v_size \n";
		exit(EXIT_FAILURE);
	}
	
	this->n_row=1;
	this->n_column=v_size;
	this->data=(double**)malloc(n_row*sizeof(double*));
	
	if(this->data==NULL){
		cout << "Matrix create: Error in data\n";
		exit(EXIT_FAILURE);
	}
	
	this->data[0]=(double*) calloc(n_column,sizeof(double));
	
	
}

Matrix::Matrix(const int n_row,const int n_column){
	if(n_row <=0 || n_column<=0){
		cout << "Matrix: Error in n_row or in n_column\n";
		exit(EXIT_FAILURE);
	}
	
	this->n_row=n_row;
	this->n_column=n_column;
	this->data=(double**)malloc(n_row*sizeof(double*));
	
	if(this->data==NULL){
		cout << "Matrix: Error in data\n";
		exit(EXIT_FAILURE);
	}
	
	for(int i=0;i<n_row;i++){
		this->data[i]=(double*)malloc(n_column*sizeof(double));
	}
	
}

double& Matrix::operator () (const int n){
	if(n<=0 || n>this->n_row*this->n_column){
		cout<<"Vector get: error in n \n";
		exit(EXIT_FAILURE);
	}
	
	return this->data[(n-1)/this->n_column][(n-1)%this->n_column];
}

double& Matrix::operator () (const int n_row,const int n_column){
	if(n_row <=0 || n_column<=0 || n_row > this->n_row || n_column > this->n_column){
	cout << "Matrix: Error in n_row or in n_column\n";
	exit(EXIT_FAILURE);
	}
	return this->data[n_row-1][n_column-1];
}

Matrix& Matrix::operator + (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sum: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + m(i,j);
		}
	}
	
	return *m_aux;
}

Matrix& Matrix::operator + (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) + s;
		}
	}
	
	return *m_aux;
}





Matrix& Matrix::operator - (Matrix &m) {
	if (this->n_row != m.n_row || this->n_column != m.n_column) {
		cout << "Matrix sub: error in n_row/n_column\n";
        exit(EXIT_FAILURE);
	}
	
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - m(i,j);
		}
	}
	
	return *m_aux;
}
Matrix& Matrix::operator - (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) - s;
		}
	}
	
	return *m_aux;
}


Matrix& Matrix::operator * (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) * s;
		}
	}
	
	return *m_aux;
}


Matrix& Matrix::operator / (double s){
	Matrix *m_aux = new Matrix(this->n_row, this->n_column);
	
    for(int i = 1; i <= this->n_row; i++) {
        for(int j = 1; j <= this->n_column; j++) {
			(*m_aux)(i,j) = (*this)(i,j) / s;
		}
	}
	
	return *m_aux;
}

ostream& operator << (ostream &o, Matrix &m) {
	for (int i = 1; i <= m.n_row; i++) {
        for (int j = 1; j <= m.n_column; j++)
			printf("%5.20lf ", m(i,j));
        o << "\n";
    }
	
    return o;
}

Matrix& zeros(const int n_row, const int n_column) {
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}
