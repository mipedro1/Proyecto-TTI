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
#include "..\include\sign_.hpp"
#include "..\include\timediff.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\IERS.hpp"
#include "../include/global.hpp"
#include "..\include\Legendre.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\gmst.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\LTC.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\gast.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"
#include "..\include\DEInteg.hpp"
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
	
	
	Matrix B(f,c);
	B(1,1) = 1; B(1,2) = 3;
	B(2,1) = 2; B(2,2) = 4;
    
    _assert(m_equals(A.transpose(), B, 1e-10));
    
    return 0;
}
int m_zeros_02() {
    int f = 2;
    int c = 2;
	
	Matrix A=zeros(3);
	
	Matrix B(1,3);
	B(1,1) = 0; B(1,2) = 0; B(1,3)=0;
    
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
int I1_sign__01() {
    
	double M=sign_(5.0,-3.0); 
    _assert( (fabs(M-(-5)))<1e-10);
    
    return 0;
}
int I1_timediff_01() {
    
	double UT1_UTC = 0.0;  
    double TAI_UTC = 37.0;       
    
    
    auto [UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
	double ut1_tai=-37;
	double utc_gps=-18;
	double ut1_gps=-18;
	double tt_utc=69.184;
	double gps_utc=18;
	
    _assert( (fabs(UT1_TAI-(ut1_tai)))<1e-10);
	_assert( (fabs(UTC_GPS-(utc_gps)))<1e-10);
	_assert( (fabs(UT1_GPS-(ut1_gps)))<1e-10);
    _assert( (fabs(TT_UTC-(tt_utc)))<1e-10);
    _assert( (fabs(GPS_UTC-(gps_utc)))<1e-10);


    
    return 0;
}
int I1_AzElPa_01() {
    
	Matrix s(3);       
    s(1,1)=1000;s(1,2)=2000;s(1,3)=500;
    
    auto [Az, El, dAds, dEds] = AzElPa(s);
	double az=0.463647609000806;
	double el=0.219987977395459;
	Matrix &dads=zeros(3);
	dads(1,1)=0.0004;dads(1,2)=-0.0002;dads(1,3)=0;
	
	Matrix &deds=zeros(3);
	deds(1,1)=-4.2591770999996e-05;deds(1,2)=-8.5183541999992e-05;deds(1,3)=0.00042591770999996;
	
    _assert( (fabs(Az-(az)))<1e-10);
	_assert( (fabs(El-(el)))<1e-10);
    _assert(m_equals(dAds,dads, 1e-10));
	_assert(m_equals(dEds,deds, 1e-10));


    
    return 0;
}

int I1_NutAngles_01() {
    
	
    double Mjd_TT=3.0;
    auto [dpsi, deps] = NutAngles (Mjd_TT);
	double DP=2.72256565175042e-05;
	double DE=3.87947551912632e-05;
    _assert( (fabs(DP-(dpsi)))<1e-10);
	_assert( (fabs(DE-(deps)))<1e-10);


    
    return 0;
}

int I1_IERS_01() {
    
	double Mjd_UTC =49746.1101504629;
    auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS (Mjd_UTC,'l');
	double x=-5.59386183152189e-07;
	double y=2.33554438440373e-06;
	double utl=0.325764698106523;
	double l=0.00272668635763815;
	double dp=-1.16881960640421e-07     ;
	double de=-2.47881680412219e-08;
	double dx=-8.41764150670523e-10     ;
	double dy=-1.56618880121342e-09;
	double tai=29;
    _assert( (fabs(x_pole-(x)))<1e-10);
	_assert( (fabs(y_pole-(y)))<1e-10);
	_assert( (fabs(utl-(UT1_UTC)))<1e-10);
	_assert( (fabs(l-(LOD)))<1e-10);
	_assert( (fabs(dp-(dpsi)))<1e-10);
	_assert( (fabs(de-(deps)))<1e-10);
	_assert( (fabs(dx-(dx_pole)))<1e-10);
	_assert( (fabs(dy-(dy_pole)))<1e-10);
	_assert( (fabs(tai-(TAI_UTC)))<1e-10);


    
    return 0;
}

int I1_Legendre_01() {
    
	
    
    auto [pnm, dpnm] = Legendre(2,2,2.5);
	Matrix &p=zeros(3,3);
	p(1,1)=1;p(1,2)=0;p(1,3)=0;
	p(2,1)=1.03658416050274;p(2,2)=-1.38762144628672;p(2,3)=0;
	p(3,1)=0.0833010473685024;p(3,2)=-1.85694887302218;p(3,3)=1.24290056661382;
	
	Matrix &d=zeros(3,3);
	d(1,1)=0;d(1,2)=0;d(1,3)=0;
	d(2,1)=-1.38762144628672;d(2,2)=-1.03658416050274;d(2,3)=0;
	d(3,1)=-3.21632979513219;d(3,2)=1.09861892024788;d(3,3)=1.85694887302218;
	
    _assert(m_equals(pnm,p, 1e-10));
	_assert(m_equals(dpnm,d, 1e-10));


    
    return 0;
}

int I1_TimeUpdate_01() {
    
	
    Matrix &P = zeros(2, 2);
	P(1,1) =1; P(1,2) =2; P(2,1)=3;P(2,2)=4;
	
	Matrix &Phi = zeros(2, 2);
	Phi(1,1) =1; Phi(1,2) =1; Phi(2,1)=1;Phi(2,2)=1;
	 
	
	Matrix &R = zeros(2,2);
	R(1,1) =10; R(1,2) =10; R(2,1)=10;R(2,2)=10;
	
	TimeUpdate(P,Phi);
	
    _assert(m_equals(P,R, 1e-10));
    
    return 0;
}

int I1_AccelHarmonic_01() {
    
	
    Matrix P(3,1);
	P(1,1) =7000e3; P(2,1) =0; P(3,1)=0;
	
	Matrix E=eye(3);
	
	 
	
	Matrix R(3,1);
	R(1,1) =-8.14576607065686; R(2,1) =-3.66267894892037e-05; R(3,1)=-5.84508413583961e-09;
	
	Matrix a=AccelHarmonic(P,E,2,2);
	
    _assert(m_equals(a,R, 1e-10));
    
    return 0;
}
int I1_Gmst_01() {
    
	
	double resultado=1.45762881716801;
	
    _assert( (fabs(gmst(60110.5)-(resultado)))<1e-10);
    
    return 0;
}

int I1_PrecMatrix_01() {
    
	
    Matrix P(3,3);
	P(1,1) =0.999988458573035; P(1,2) =-0.0044064693314773; P(1,3)=-0.00191461451888747;
	P(2,1) =0.00440646933120376; P(2,2) =0.999990291457991; P(2,3)=-4.21851229096372e-06;
	P(3,1) =0.00191461451951702; P(3,2) =-4.21822655524287e-06; P(3,3)=0.999998167115044;
	
	Matrix R=PrecMatrix(51544.0,58741.0);
	
    _assert(m_equals(P,R, 1e-10));
    
    return 0;
}
int I1_PoleMatrix_01() {
    
	
    Matrix P(3,3);
	P(1,1) =0.999998000000667; P(1,2) =1.99999833333384e-06; P(1,3)=0.00199999766666768;
	P(2,1) =0; P(2,2) =0.999999500000042; P(2,3)=-0.000999999833333342;
	P(3,1) =-0.00199999866666693; P(3,2) =0.000999997833334342; P(3,3)=0.999997500001708;
	
	Matrix R=PoleMatrix(0.002,0.001);
	
    _assert(m_equals(P,R, 1e-10));
    
    return 0;
}

int I1_NutMatrix_01() {
    
	
    Matrix P(3,3);
	P(1,1) =0.999999997583483; P(1,2) =6.37847604742929e-05; P(1,3)=2.76502946008031e-05;
	P(2,1) =-6.37842142060874e-05; P(2,2) =0.999999997770622; P(2,3)=-1.97567566176948e-05;
	P(3,1) =-2.76515547191487e-05; P(3,2) =1.97549929176755e-05; P(3,3)=0.999999999422566;
	
	Matrix R=NutMatrix(59580.0);
	
    _assert(m_equals(P,R, 1e-10));
    
    return 0;
}

int I1_LTC_01() {
    
	
    Matrix P(3,3);
	P(1,1) =0.961277227548362; P(1,2) =0.275583184895113; P(1,3)=0;
	P(2,1) =-0.179676221403558; P(2,2) =0.626738746897472; P(2,3)=0.758231494069936;
	P(3,1) =0.208955850023573; P(3,2) =-0.7288706684594; P(3,3)=0.651985430359048;
	
	Matrix R=LTC(-1.2916,0.7102);
	
    _assert(m_equals(P,R, 1e-10));
    
    return 0;
}
int I1_EqnEquinox_01() {
    
	
	double resultado=-6.19313954901517e-05;
	
    _assert( (fabs(EqnEquinox(51544.0)-(resultado)))<1e-10);
    
    return 0;
}
int I1_Gast_01() {
    
	
	double resultado=1.74470523193512;
	
    _assert( (fabs(gast(51544.0)-(resultado)))<1e-10);
    
    return 0;
}
int I1_Jpl_Ep_01() {
    
	
    
    auto [r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun] = JPL_Eph_DE430(59580.0);
	Matrix r_m(3,1);
	r_m(1,1)=79837647856.2735;
	r_m(2,1)=-136261500272.064;
	r_m(3,1)=-64981938909.0024;
	
	Matrix r_v(3,1);
	r_v(1,1)=15975346745.7943;
	r_v(2,1)=-35365322602.8999;
	r_v(3,1)=-13084365467.2929;
	
	Matrix r_e(3,1);
	r_e(1,1)=-27411369555.7101;
	r_e(2,1)=133274848816.261;
	r_e(3,1)=57802495046.4152;
	
	Matrix r_mar(3,1);
	r_mar(1,1)=-103539704034.338;
	r_mar(2,1)=-306659028736.292;
	r_mar(3,1)=-133814565360.54;
	
	Matrix r_j(3,1);
	r_j(1,1)=722966922027.3;
	r_j(2,1)=-373377540853.151;
	r_j(3,1)=-177648967214.498;
	
	Matrix r_s(3,1);
	r_s(1,1)=1067340429956.83;
	r_s(2,1)=-1093570901121.93;
	r_s(3,1)=-499244305007.512;
	
	Matrix r_u(3,1);
	r_u(1,1)=2179981852798.6;
	r_u(2,1)=1725296913450.38;
	r_u(3,1)=725756527252.905;
	
	Matrix r_nep(3,1);
	r_nep(1,1)=4459200999344.89;
	r_nep(2,1)=-658649518957.463;
	r_nep(3,1)=-383177121895.469;
	
	Matrix r_p(3,1);
	r_p(1,1)=2289056144424.15;
	r_p(2,1)=-4312612949408.77;
	r_p(3,1)=-2043473055509.64;
	
	Matrix r_moo(3,1);
	r_moo(1,1)=-91868548.3240412;
	r_moo(2,1)=-315040557.002235;
	r_moo(3,1)=-145304426.66531;
	
	Matrix r_su(3,1);
	r_su(1,1)=26127801143.8142;
	r_su(2,1)=-132825709286.88;
	r_su(3,1)=-57579560410.2985;
	
    _assert(m_equals(r_Mercury,r_m, 1e-2));
	_assert(m_equals(r_Venus,r_v, 1e-2));
	_assert(m_equals(r_Earth,r_e, 1e-2));
	_assert(m_equals(r_Mars,r_mar, 1e-2));
	_assert(m_equals(r_Jupiter,r_j, 1e-2));
	_assert(m_equals(r_Saturn,r_s, 1e-2));
	_assert(m_equals(r_Uranus,r_u, 1e-2));
	_assert(m_equals(r_Neptune,r_nep, 1e-2));
	_assert(m_equals(r_Pluto,r_p, 1e-2));
	_assert(m_equals(r_Moon,r_moo, 1e-2));
	_assert(m_equals(r_Sun,r_su, 1e-2));


    
    return 0;
}
int I1_MeasUpdate_01() {
    
	Matrix &x=zeros(2,2);
	x(1,1)=1;x(1,2)=2;
	x(2,1)=1;x(2,2)=2;
    Matrix &G=eye(2);
	Matrix &P=zeros(2,2);
	P(1,1)=0.5; P(1,2)=0;
	P(2,1)=0; P(2,2)=0.5;
	
    auto [K, xx, pp] = MeasUpdate(x,3,2.5,0.2,G,P,2);
	Matrix &kk=zeros(2,2);
	kk(1,1)=0.931034482758621;kk(1,2)=-0.0689655172413793;
	kk(2,1)=-0.0689655172413793;kk(2,2)=0.931034482758621;
	
	Matrix &x_result=zeros(2,2);
	x_result(1,1)=1.46551724137931;x_result(1,2)=1.96551724137931;
	x_result(2,1)=0.96551724137931;x_result(2,2)=2.46551724137931;
	
	Matrix &p_result=zeros(2,2);
	p_result(1,1)=0.0344827586206897;p_result(1,2)=0.0344827586206897;
	p_result(2,1)=0.0344827586206897;p_result(2,2)=0.0344827586206897;
	
    _assert(m_equals(K,kk, 1e-10));
	_assert(m_equals(x_result,xx, 1e-10));
	_assert(m_equals(p_result,pp, 1e-10));


    
    return 0;
}
int I1_G_Accel_01() {
    
	
	Matrix &x=zeros(3,3);
	x(1,1)=2.33052264420053e-06;x(1,2)=2.09290362818138e-11;x(1,3)=3.5527136788005e-15;
	x(2,1)=2.0929593988635e-11;x(2,2)=-1.16369909716371e-06;x(2,3)=5.47310677681545e-15;
	x(3,1)=3.34004807862169e-15;x(3,2)=5.47310880258077e-15;x(3,3)=-1.16682354671017e-06;
	
	
	Matrix P(3,1);
	P(1,1) =7000e3; P(2,1) =0; P(3,1)=0;
	
	Matrix &E=eye(3);
	
	
	
    _assert(m_equals(G_AccelHarmonic(P,E,2,2),x, 1e-10));


    
    return 0;
}
int I1_GHAMatrix_01() {
    
	
	Matrix &x=zeros(3,3);
	x(1,1)=-0.184404138193381;x(1,2)=0.982850504307322;x(1,3)=0;
	x(2,1)=-0.982850504307322;x(2,2)=-0.184404138193381 ;x(2,3)=0;
	x(3,1)=0;x(3,2)=0;x(3,3)=1;
	
	
	
	
	
    _assert(m_equals(GHAMatrix(59580),x, 1e-10));


    
    return 0;
}

int I1_Accel_01() {
    
	
	Matrix &x=zeros(6,1);
	x(1,1)=0;
	x(2,1)=7.5;
	x(3,1)=0;
	x(4,1)=-1.43703755432425e+60;
	x(5,1)=-2.25906002132808e+59;
	x(6,1)=-3.71905163213436e+60;
	
	
	Matrix &p=zeros(6,1);
	p(1,1)=7000;
	p(2,1)=0;
	p(3,1)=0;
	p(4,1)=0;
	p(5,1)=7.5;
	p(6,1)=0;
    _assert(m_equals(Accel(51544.0,p),x, 1e+51));


    
    return 0;
}

int I1_VarEqn_01() {
    
	
	Matrix &x=zeros(42,1);
	x(1,1)=5394.06842166351;
	x(2,1)=-2365.21337882342;
	x(3,1)=-7061.84554200295;
	x(4,1)=-5.1348367854085;
	x(5,1)=-2.97717622353621;
	x(6,1)=-3.70591776714203;
	x(10,1)=5.70032034019619e-07;
	x(11,1)=8.67651591462959e-07;
	x(12,1)=1.08169353829624e-06;
	x(16,1)=8.67651589686602e-07;
	x(17,1)=-4.23359106882515e-07;
	x(18,1)=6.27183701418232e-07;
	x(22,1)=1.08169353918441e-06;
	x(23,1)=6.27183704526857e-07;
	x(24,1)=-1.46672927137104e-07;
	x(25,1)=1;
	x(32,1)=1;
	x(39,1)=1;
	
	
	Matrix &p=zeros(42,1);
	p(1,1)=5542555.93722861;
	p(2,1)=3213514.8673492;
	p(3,1)=3990892.97587685;
	p(4,1)=5394.06842166351;
	p(5,1)=-2365.21337882342;
	p(6,1)=-7061.84554200295;
	p(7,1)=1;
	p(14,1)=1;
	p(21,1)=1;
	p(28,1)=1;
	p(35,1)=1;
	p(42,1)=1;
    _assert(m_equals(VarEqn(51544.0,p),x, 1e+51));


    
    return 0;
}
int I1_DeInteg_01() {
    
	
	Matrix &x=zeros(6,1);
	x(1,1)=5542555.89427451;
	x(2,1)=3213514.83814162;
	x(3,1)=3990892.92789074;
	x(4,1)=5394.06894044389;
	x(5,1)=-2365.21290574021;
	x(6,1)=-7061.8448137347;
	
	Matrix &Y0_apr=zeros(6,1);
	Y0_apr(1,1)=6221397.62857869;
	Y0_apr(2,1)=2867713.77965738;
	Y0_apr(3,1)=3006155.98509949;
	Y0_apr(4,1)=4645.04725161806;
	Y0_apr(5,1)=-2752.21591588204;
	Y0_apr(6,1)=-7507.99940987031;
	
	double t_aux=0.0;
	double obs=-134.999991953373;
	Y0_apr=DEInteg(Accel,t_aux,obs,1e-13,1e-6,6,Y0_apr);
	cout<<"Y0_apr\n"<<Y0_apr<<endl;
    _assert(m_equals(Y0_apr,x, 1e-10));


    
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
	_verify(I1_sign__01);
	_verify(I1_timediff_01);
	_verify(I1_AzElPa_01);
	_verify(I1_IERS_01);
	_verify(I1_NutAngles_01);
	_verify(I1_Legendre_01);
	_verify(I1_TimeUpdate_01);
	_verify(I1_Gmst_01);
	_verify(I1_PrecMatrix_01);
	_verify(I1_AccelHarmonic_01);
	_verify(I1_PoleMatrix_01);
	_verify(I1_NutMatrix_01);
	_verify(I1_LTC_01);
	_verify(I1_EqnEquinox_01);
	_verify(I1_Gast_01);
	_verify(I1_Jpl_Ep_01);
	_verify(I1_MeasUpdate_01);
	_verify(I1_G_Accel_01);
	_verify(I1_GHAMatrix_01);
	_verify(I1_Accel_01);
	_verify(I1_VarEqn_01);
	_verify(I1_DeInteg_01);
	
	

    return 0;
}


int main()
{
	eop19620101(21413);
	GGM03S(181);
	DE430Coeff(2285,1020);
	AuxParam.Mjd_UTC=49746.1112847221;
	AuxParam.Mjd_TT=49746.1170623147;
	AuxParam.n      = 20;
	AuxParam.m      = 20;
	AuxParam.sun     = 1;
	AuxParam.moon    = 1;
	AuxParam.planets = 1;
    int result = all_tests();

    if (result == 0)
        printf("PASSED\n");

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}
