#include "../include/matrix.hpp"


Matrix::Matrix(){
	this->n_row=0;
	this->n_column=0;
	this->data= nullptr;
}

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
		cout << "Matrix\n"<<*this<<endl;
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
	cout << "Matrix\n"<<*this<<endl;
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


Matrix& Matrix::operator * (Matrix &m) {
    
	if (this->n_column != m.n_row) {
        cout << "Matrix multiplication: error in dimensions\n";
        exit(EXIT_FAILURE);
    }

    Matrix *m_aux = new Matrix(this->n_row, m.n_column);

    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= m.n_column; j++) {
            (*m_aux)(i, j) = 0;
            for (int k = 1; k <= this->n_column; k++) {
                (*m_aux)(i, j) += (*this)(i, k) * m(k, j); 
            }
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


Matrix& Matrix::operator / (Matrix &m) {
    if (this->n_column != m.n_row) {
        cout << "Matrix division: error in dimensions\n";
        exit(EXIT_FAILURE);
    }

    if (m.n_row != m.n_column) {
        cout << "Matrix division: Matrix B must be square to compute its inverse.\n";
        exit(EXIT_FAILURE);
    }

    if (m.det() == 0) {
        cout << "Matrix division: Matrix B is singular and does not have an inverse.\n";
        exit(EXIT_FAILURE);
    }

    Matrix inversa = m.inv();  

    
    Matrix *m_aux = new Matrix(this->n_row, this->n_column);

    *m_aux = (*this) * inversa;


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

Matrix& Matrix::operator = (Matrix& m) {
    if (this == &m) {
        return *this;
    }

    //liberar memoria
    for (int i = 0; i < this->n_row; i++) {
        delete[] this->data[i];
    }
    delete[] this->data; 

    
    this->n_row = m.n_row;
    this->n_column = m.n_column;
    this->data = new double*[this->n_row];

    for (int i = 0; i < this->n_row; ++i) {
        this->data[i] = new double[this->n_column];
    }

    
    for (int i = 0; i < this->n_row; ++i) {
        for (int j = 0; j < this->n_column; ++j) {
            this->data[i][j] = m.data[i][j];
        }
    }

    
    return *this;
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
	if (n_row <= 0 || n_column<=0) {
        cout << "Error in zeros: n_row and n_column must be positive\n";
        exit(EXIT_FAILURE);
    }
	Matrix *m_aux = new Matrix(n_row, n_column);
	
	for(int i = 1; i <= n_row; i++) {
		for(int j = 1; j <= n_column; j++) {
			(*m_aux)(i,j) = 0;
		}
	}
	
	return (*m_aux);
}

Matrix& Matrix::inv() {
	if (this->n_row != this->n_column) {
		cout << "Matrix inv: matrix must be square\n";
        exit(EXIT_FAILURE);
	}

    int n = this->n_row;
    
    // Crear la matriz aumentada [A | I]
    Matrix* aumentada = new Matrix(n, 2 * n);

    // Llenar la parte izquierda con la matriz original y la derecha con la identidad
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*aumentada)(i, j) = (*this)(i, j);        // Copiar A
            (*aumentada)(i, j + n) = (i == j) ? 1 : 0; // Copiar I
        }
    }

    // Gauss-Jordan
    for (int i = 1; i <= n; i++) {
        // Encontrar el pivote (elemento no nulo más grande en la columna)
        double pivot = (*aumentada)(i, i);
        if (pivot == 0) {
            cout << "The matrix is singular and does not have an inversa\n";
            exit(EXIT_FAILURE);
        }

        // Dividir la fila por el pivote
        for (int j = 1; j <= 2 * n; j++) {
            (*aumentada)(i, j) /= pivot;
        }

        // Hacer ceros en la columna del pivote
        for (int j = 1; j <= n; j++) {
            if (i != j) {
                double factor = (*aumentada)(j, i);
                for (int k = 1; k <= 2 * n; k++) {
                    (*aumentada)(j, k) -= factor * (*aumentada)(i, k);
                }
            }
        }
    }

    // Extraer la parte derecha como la matriz inversa
    Matrix* inversa = new Matrix(n, n);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= n; j++) {
            (*inversa)(i, j) = (*aumentada)(i, j + n); // Parte derecha de la matriz aumentada
        }
    }

    return *inversa;
}


double Matrix::det() const {
    if (this->n_row != this->n_column) {
		cout << "Matrix inv: matrix must be square\n";
        exit(EXIT_FAILURE);
	}

    // matrices 2x2
    if (n_row == 2 && n_column == 2) {
        return data[0][0] * data[1][1] - data[0][1] * data[1][0];
    }

    // matrices 3x3
    if (n_row == 3 && n_column == 3) {
        return data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1])
             - data[0][1] * (data[1][0] * data[2][2] - data[1][2] * data[2][0])
             + data[0][2] * (data[1][0] * data[2][1] - data[1][1] * data[2][0]);
    }

    // matrices más grandes (implementación recursiva de determinante)
    if (n_row > 3) {
        double det_value = 0;
        for (int i = 0; i < n_column; i++) {
            Matrix submatrix(n_row - 1, n_column - 1);
            // Crear submatriz eliminando la primera fila y la columna i
            for (int j = 1; j < n_row; j++) {
                int sub_col = 0;
                for (int k = 0; k < n_column; k++) {
                    if (k != i) {
                        submatrix(j - 1, sub_col++) = data[j][k];
                    }
                }
            }
            // Cofactor y determinante recursivo
            det_value += (i % 2 == 0 ? 1 : -1) * data[0][i] * submatrix.det();
        }
        return det_value;
    }
	return 0;
}
Matrix& eye(int n) {
    if (n <= 0) {
        cout << "Matrix eye: Invalid size n\n";
        exit(EXIT_FAILURE);
    }
	
	Matrix *m_aux = new Matrix(n, n);
	
	for(int i = 1; i <= n; i++) {
		for(int j = 1; j <= n; j++) {
			if (i == j) {
                (*m_aux)(i,j) = 1;
            } else {
                (*m_aux)(i,j) = 0;
            }
			
		}
	}
	
	return (*m_aux);
}

Matrix& Matrix::transpose() {
    Matrix* m_aux = new Matrix(this->n_column, this->n_row);

    
    for (int i = 1; i <= this->n_row; i++) {
        for (int j = 1; j <= this->n_column; j++) {
            (*m_aux)(j, i) = (*this)(i, j);  
        }
    }

    

    
    return (*m_aux);
}

Matrix& zeros(int n) {
	if (n <= 0) {
        cout << "Error in zeros: n must be positive\n";
        exit(EXIT_FAILURE);
    }
    Matrix *m_aux = new Matrix(1, n);  
    
    
    for(int j = 1; j <= n; j++) {
        (*m_aux)(1, j) = 0; 
    }
    
    return (*m_aux);
}

double norm(Matrix& m) {
    double sum = 0.0;

    for (int i = 1; i <= m.n_row; ++i) {
        for (int j = 1; j <= m.n_column; ++j) {
            sum += pow(m(i, j), 2);
        }
    }

    return sqrt(sum);  
}

double Matrix::dot(Matrix& m) {
    if (this->n_row != 1 || m.n_row != 1 || this->n_column != m.n_column) {
        cout << "Dot product: Vectors must have the same number of columns and must be row vectors.\n";
        exit(EXIT_FAILURE);
    }

    double result = 0.0;

    for (int j = 1; j <= this->n_column; j++) {
        result += (*this)(1, j) * m(1, j); 		
    }

    return result;
}

Matrix& Matrix::cross(Matrix& m){
    if (this->n_row != 1 || this->n_column != 3 || m.n_row != 1 || m.n_column != 3) {
        cout << "Cross product: Both matrices must be row vectors of dimension 3.\n";
        exit(EXIT_FAILURE);
    }

    
    Matrix* result = new Matrix(1, 3);

    
    (*result)(1, 1) = (*this)(1, 2) * m(1, 3) - (*this)(1, 3) * m(1, 2);  
    (*result)(1, 2) = (*this)(1, 3) * m(1, 1) - (*this)(1, 1) * m(1, 3);  
    (*result)(1, 3) = (*this)(1, 1) * m(1, 2) - (*this)(1, 2) * m(1, 1);  

    return *result;
	
}

Matrix& Matrix::extract_vector(int indiceInicio, int indiceFinal) {
    if (indiceInicio < 1 || indiceFinal > this->n_column || indiceInicio > indiceFinal) {
        cout << "Error: Invalid index range\n";
        exit(EXIT_FAILURE);
    }

    int tamañoVector = indiceFinal - indiceInicio + 1;
    Matrix* subvector = new Matrix(1, tamañoVector);

    for (int i = indiceInicio; i <= indiceFinal; i++) {
        (*subvector)(1, i - indiceInicio + 1) = (*this)(1, i);
    }

    return *subvector;
}

Matrix& Matrix::union_vector(Matrix& m) {
    if (this->n_row != 1 || m.n_row != 1) {
		 cout << "Error: Invalid index range\n";
        exit(EXIT_FAILURE);
    }

    int n_column = this->n_column + m.n_column;
    Matrix* result = new Matrix(1, n_column);

    for (int i = 0; i < this->n_column; i++) {
        (*result)(1, i + 1) = (*this)(1, i + 1);
    }

    for (int i = 0; i < m.n_column; i++) {
        (*result)(1, this->n_column + i + 1) = m(1, i + 1);
    }

    return *result;
}
Matrix& Matrix::extract_row(int n) {
    if (n <= 0 || n > this->n_row) {
        exit(EXIT_FAILURE);
    }

    Matrix* result = new Matrix(1, this->n_column);

    for (int i = 0; i < this->n_column; i++) {
        (*result)(1, i + 1) = (*this)(n, i + 1);
    }

    return *result;
}
Matrix& Matrix::extract_column(int n) {
    if (n <= 0 || n > this->n_column) {
        exit(EXIT_FAILURE);
    }

    Matrix* result = new Matrix(1, this->n_row);

    for (int i = 0; i < this->n_row; i++) {
        (*result)(1, i + 1) = (*this)(i + 1, n);
    }

    return *result;
}

Matrix& Matrix::assign_row(int n, Matrix& m) {
    if (n <= 0 || n > this->n_row || m.n_column != this->n_column) {
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < this->n_column; i++) {
        (*this)(n, i + 1) = m(1, i + 1);
    }

    return *this;
}

Matrix& Matrix::assign_column(int n, Matrix& m) {
    if (n <= 0 || n > this->n_column || m.n_row != this->n_row) {
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < this->n_row; i++) {
        (*this)(i + 1, n) = m(i + 1, 1);
    }

    return *this;
}


