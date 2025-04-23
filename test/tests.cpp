#include "..\include\matrix.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Cheb3D.hpp"
#include "..\include\EccAnom.hpp"
#include "..\include\Frac.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\Mjday.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\Position.hpp"
#include <cstdio>
#include <cmath>

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

int m_equals(Matrix A, Matrix B, double p) {
	if (A.n_row != B.n_row || A.n_column != B.n_column)
		return 0;
	else
		for(int i = 1; i <= A.n_row; i++)
			for(int j = 1; j <= A.n_column; j++)
				if(fabs(A(i,j)-B(i,j)) > p) {
					printf("%2.20lf %2.20lf\n",A(i,j),B(i,j));
					return 0;
				}
	
	return 1;
}

int m_sum_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = 8; C(2,2) = -3; C(2,3) = 1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = -2; C(3,3) = 0; C(3,4) = 7;
	
	Matrix R = A + B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}
int m_sum_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	double suma=3;
	
	Matrix C(f, c);
	C(1,1) = 3; C(1,2) =  5; C(1,3) = 11; C(1,4) = 3;
	C(2,1) = 4; C(2,2) = 2; C(2,3) = 3; C(2,4) = 3;
	C(3,1) = 3; C(3,2) = 4; C(3,3) = 3; C(3,4) = 8;
	
	Matrix R = A + suma;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}
int m_res_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = -2; C(1,2) =  2; C(1,3) = 8; C(1,4) = 0;
	C(2,1) = -6; C(2,2) = 1; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 4; C(3,3) = 0; C(3,4) = 3;
	
	Matrix R = A - B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}
int m_res_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	double num=3;
	
	Matrix C(f, c);
	C(1,1) = -3; C(1,2) =  -1; C(1,3) = 5; C(1,4) = -3;
	C(2,1) = -2; C(2,2) = -4; C(2,3) = -3; C(2,4) = -3;
	C(3,1) = -3; C(3,2) = -2; C(3,3) = -3; C(3,4) = 2;
	
	Matrix R = A - num;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}
int m_mul_01() {
    int f = 4;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	A(4,1) = 0; A(4,2) =  1; A(4,3) = 0; A(4,4) = 5;
	
	Matrix B(f, c);
	B(1,1) = 2; B(1,2) =  0; B(1,3) = 0; B(1,4) = 0;
	B(2,1) = 7; B(2,2) = -2; B(2,3) = 1; B(2,4) = 0;
	B(3,1) = 0; B(3,2) = -3; B(3,3) = 0; B(3,4) = 2;
	B(4,1) = 0; B(4,2) = -3; B(4,3) = 0; B(4,4) = 2;
	
	Matrix C(f, c);
	C(1,1) = 14; C(1,2) =  -28; C(1,3) = 2; C(1,4) = 16;
	C(2,1) = -5; C(2,2) = 2; C(2,3) = -1; C(2,4) = 0;
	C(3,1) = 7; C(3,2) = -17; C(3,3) = 1; C(3,4) = 10;
	C(4,1) = 7; C(4,2) = -17; C(4,3) = 1; C(4,4) = 10;
	
	Matrix R = A * B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}
int m_mul_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	double num=2;
	
	Matrix C(f, c);
	C(1,1) = 0; C(1,2) =  4; C(1,3) = 16; C(1,4) = 0;
	C(2,1) = 2; C(2,2) = -2; C(2,3) = 0; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 2; C(3,3) = 0; C(3,4) = 10;
	
	Matrix R = A * num;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}
int m_div_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) =  2;
	A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(f, c);
	B(1,1) = 1; B(1,2) =  1;
	B(2,1) = 1; B(2,2) = 2; 
	
	Matrix C(f, c);
	C(1,1) = 0; C(1,2) = 1;
	C(2,1) = 2; C(2,2) = 1; 
	
	Matrix R = A / B;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_div_02() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) =  2; A(1,3) = 8; A(1,4) = 0;
	A(2,1) = 1; A(2,2) = -1; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) =  1; A(3,3) = 0; A(3,4) = 5;
	
	double num=2;
	
	Matrix C(f, c);
	C(1,1) = 0; C(1,2) =  1; C(1,3) = 4; C(1,4) = 0;
	C(2,1) = 0.5; C(2,2) = -0.5; C(2,3) = 0; C(2,4) = 0;
	C(3,1) = 0; C(3,2) = 0.5; C(3,3) = 0; C(3,4) = 2.5;
	
	Matrix R = A / num;
    
    _assert(m_equals(C, R, 1e-10));
    
    return 0;
}

int m_zeros_01() {
    int f = 3;
    int c = 4;
	
	Matrix A(f, c);
	A(1,1) = 0; A(1,2) = 0; A(1,3) = 0; A(1,4) = 0;
	A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0;
	A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 0;
	
	Matrix B = zeros(3, 4);
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}
int m_inv_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 1;
	A(2,1) = 1; A(2,2) = 2;
	
	Matrix inv=A.inv();
	
	Matrix B(f, c);
    B(1,1) = 2; B(1,2) = -1;
	B(2,1) = -1; B(2,2) = 1;
    _assert(m_equals( inv, B, 1e-10));
    
    return 0;
}

int m_opigual_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) =  2;
	A(2,1) = 3; A(2,2) = 4;
	
	Matrix B = A;
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}
int m_det_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 3; A(2,2) = 4;
	
	double det=-2;
    
    _assert((A.det()-det)<1e-10);
    
    return 0;
}

int m_eye_01() {
    int f = 2;
    int c = 2;
	
	Matrix A=eye(2);
	
	Matrix B(f,c);
	B(1,1) = 1; B(1,2) = 0;
	B(2,1) = 0; B(2,2) = 1;
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}

int m_trans_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 3; A(2,2) = 4;
	
	A.transpose();
	Matrix B(f,c);
	B(1,1) = 1; B(1,2) = 3;
	B(2,1) = 2; B(2,2) = 4;
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}
int m_zeros_02() {
    int f = 2;
    int c = 2;
	
	Matrix A=zeros(3);
	
	Matrix B(3,3);
	B(1,1) = 0; B(1,2) = 0; B(1,3)=0;
	B(2,1) = 0; B(2,2) = 0; B(2,3)=0;
	B(3,1) = 0; B(3,2) = 0; B(3,3)=0;
    
    _assert(m_equals(A, B, 1e-10));
    
    return 0;
}
int m_norm_01() {
    int f = 1;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = 3;
	
	
	
    
    _assert((norm(A)-3.74165738677394)<1e-10);
    
    return 0;
}
int m_dot_01() {
    int f = 1;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; 
	
	Matrix B(f, c);
	B(1,1) = 3; B(1,2) = 2;
	
    
    _assert((A.dot(B)-7)<1e-10);
    
    return 0;
}
int m_cross_01() {
    int f = 1;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = -4;
	Matrix B(f, c);
	B(1,1) = 1; B(1,2) = 2; B(1,3) = 3;
	
	Matrix C(f, c);
	C(1,1) = 14; C(1,2) = -7; C(1,3) = 0;
	
	
    
    _assert(m_equals(C, A.cross(B), 1e-10));
    
    return 0;
}
int m_extract_01() {
    int f = 1;
    int c = 3;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; A(1,3) = -4;
	
	
	Matrix C(f, 2);
	C(1,1) = 2; C(1,2) = -4; 
	
	
    
    _assert(m_equals(C, A.extract_vector(2,3), 1e-10));
    
    return 0;
}
int m_union_01() {
    int f = 1;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2; 
	
	Matrix B(f, c);
	B(1,1) = 1; B(1,2) = 2; 
	
	
	Matrix C(f, 4);
	C(1,1) = 1; C(1,2) = 2;C(1,3) = 1; C(1,4) = 2;  
	
	
    
    _assert(m_equals(C, A.union_vector(B), 1e-10));
    
    return 0;
}
int m_exrow_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 3; A(2,2) = 4;
	
	
	
	
	Matrix C(1,2);
	C(1,1) = 3; C(1,2) = 4;
	
	
    
    _assert(m_equals(C, A.extract_row(2), 1e-10));
    
    return 0;
}
int m_excol_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 3; A(2,2) = 4;
	
	
	
	
	Matrix C(1,2);
	C(1,1) = 2; C(1,2) = 4;
	
	
    
    _assert(m_equals(C, A.extract_column(2), 1e-10));
    
    return 0;
}
int m_assrow_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(1,2);
	B(1,1)= 2; B(1,2)= 2;
	
	
	
	Matrix C(f,c);
	C(1,1) = 1; C(1,2) = 2;
	C(2,1) = 2; C(2,2) = 2;
	
	
    
    _assert(m_equals(C, A.assign_row(2,B), 1e-10));
    
    return 0;
}
int m_asscol_01() {
    int f = 2;
    int c = 2;
	
	Matrix A(f, c);
	A(1,1) = 1; A(1,2) = 2;
	A(2,1) = 3; A(2,2) = 4;
	
	Matrix B(2,1);
	B(1,1)= 2;
	B(2,1)= 2;
	
	
	
	Matrix C(f,c);
	C(1,1) = 1; C(1,2) = 2;
	C(2,1) = 3; C(2,2) = 2;
	
	
    
    _assert(m_equals(C, A.assign_column(2,B), 1e-10));
    
    return 0;
}
int I1_R_x_01() {
    
	
	
	Matrix A(3, 3);
	A(1,1) = 1; A(1,2) = 0; A(1,3)=0;
	A(2,1) = 0; A(2,2) = -0.989992496600445 ; A(2,3)=0.141120008059867;
	A(3,1) = 0; A(3,2) = -0.141120008059867; A(3,3)= -0.989992496600445;
	        
	       
	Matrix B= R_x(3);
	
	
    _assert(m_equals(A,B, 1e-10));
    
    return 0;
}
int I1_R_y_01() {
    
	
	
	Matrix A(3, 3);
	A(1,1) = -0.989992496600445; A(1,2) = 0; A(1,3)=-0.141120008059867;
	A(2,1) = 0; A(2,2) =  1; A(2,3)=0;
	A(3,1) = 0.141120008059867; A(3,2) =0 ; A(3,3)= -0.989992496600445;
	        
	       
	Matrix B= R_y(3);
	
	
    _assert(m_equals(A,B, 1e-10));
    
    return 0;
}
int I1_R_z_01() {
    
	
	
	Matrix A(3, 3);
	A(1,1) = -0.989992496600445; A(1,2) = 0.141120008059867; A(1,3)=0;
	A(2,1) = -0.141120008059867; A(2,2) = -0.989992496600445 ; A(2,3)=0;
	A(3,1) = 0; A(3,2) = 0; A(3,3)= 1;
	        
	       
	Matrix B= R_z(3);
	
	
    _assert(m_equals(A,B, 1e-10));
    
    return 0;
}
int I1_AccelPointMass_01() {
    
	
	
	Matrix A(1, 3);
	A(1,1) =1; A(1,2) =2; A(1,3)=3;
	
	Matrix B(1, 3);
	B(1,1) =4; B(1,2) =5; B(1,3)=6;
	 
	double g=9.876;
	       
	Matrix C=AccelPointMass(A, B,g);
	
	Matrix D(1,3);
	D(1,1)=0.152715682717329; D(1,2)=0.138099128780766; D(1,3)=0.123482574844203;
	
	
    _assert(m_equals(C,D, 1e-10));
    
    return 0;
}
int I1_Cheb3D_01() {
    
	
	
	Matrix A(1, 5);
	A(1,1) =0.0; A(1,2) =1.0; A(1,3)=0.5;A(1,4)=0.25;A(1,5)=0.1;
	
	Matrix B(1, 5);
	B(1,1) =0.0; B(1,2) =0.9; B(1,3)=0.4;B(1,4)=0.3;B(1,5)=0.2;
	
	Matrix C(1, 5);
	C(1,1) =0.0; C(1,2) =0.7; C(1,3)=0.6;C(1,4)=0.1;C(1,5)=0.05;
	 
	
	       
	Matrix E=Cheb3D(0.5,5,0.0,1.0,A,B,C);
	
	Matrix D(1,3);
	D(1,1)=-0.4000; D(1,2)=-0.2000; D(1,3)= -0.5500;
	
	
    _assert(m_equals(E,D, 1e-10));
    
    return 0;
}

int I1_EccAnom_01() {
    
	       
	double E=EccAnom(0.5,0);
	
	
	
	
    _assert( (fabs(E-0.5000))<1e-10);
    
    return 0;
}
int I1_Frac_01() {
    
	       
	double F=Frac(4.596);
	
	
	
	
    _assert( (fabs(F-0.596))<1e-10);
    
    return 0;
}
int I1_MeanObl_01() {
    
	double M=MeanObliquity(58000.0); 
    _assert( (fabs(M-0.40905268985035))<1e-10);
    
    return 0;
}

int I1_Mjday_01() {
    
	double M=Mjday(2000, 1, 1, 0, 0, 0.0); 
    _assert( (fabs(M-51544))<1e-10);
    
    return 0;
}

int I1_MjdayTDB_01() {
    
	double M=Mjday_TDB(59900.0); 
    _assert( (fabs(M-59899.9999999855))<1e-10);
    
    return 0;
}
int I1_Position_01() {
    
	Matrix M(3);
	M(1,1) =6378136.3; M(1,2) =0; M(1,3)=0;
	Matrix r= Position(0,0,0);
    _assert(m_equals(M,r, 1e-10));
    
    return 0;
}

int all_tests()
{
    _verify(m_sum_01);
	_verify(m_sum_02);
	_verify(m_res_01);
	_verify(m_res_02);
	_verify(m_mul_01);
	_verify(m_mul_02);
	_verify(m_div_01);
	_verify(m_div_02);
    _verify(m_zeros_01);
	_verify(m_inv_01);
	_verify(m_opigual_01);
	_verify(m_det_01);
	_verify(m_eye_01);
	_verify(m_trans_01);
	_verify(m_zeros_02);
	_verify(m_norm_01);
	_verify(m_dot_01);
	_verify(m_cross_01);
	_verify(m_extract_01);
	_verify(m_union_01);
	_verify(m_exrow_01);
	_verify(m_excol_01);
	_verify(m_assrow_01);
	_verify(m_asscol_01);
	_verify(I1_R_x_01);
	_verify(I1_R_y_01);
	_verify(I1_R_z_01);
	_verify(I1_AccelPointMass_01);
	_verify(I1_Cheb3D_01);
	_verify(I1_EccAnom_01);
	_verify(I1_Frac_01);
	_verify(I1_MeanObl_01);
	_verify(I1_Mjday_01);
	_verify(I1_MjdayTDB_01);
	_verify(I1_Position_01);

    return 0;
}


int main()
{
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
