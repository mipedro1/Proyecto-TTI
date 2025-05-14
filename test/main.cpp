#include "../include/matrix.hpp"
#include "../include/R_x.hpp"
#include "../include/global.hpp"

int main(){
	eop19620101(4); // c=21413
	cout<< "eopdata\n"<<eopdata<<"\n";
	
	GGM03S(3); // n=181
	cout<< "Cnm\n"<<Cnm<<"\n";
	cout<< "Snm\n"<<Snm<<"\n";
	
	DE430Coeff(10,5); // f=2285 c=1020
	cout<< "PC\n"<<PC<<"\n";
	
	AuxParam.Mjd_UTC=49746.1112847221;
	AuxParam.Mjd_TT=49746.1170623147;
	AuxParam.n      = 20;
	AuxParam.m      = 20;
	AuxParam.sun     = 1;
	AuxParam.moon    = 1;
	AuxParam.planets = 1;
	
	Matrix aux= R_x(3);
	cout << "aux\n" << aux << "\n";
	
	Matrix v(3);
	v(2)=5;
	cout<<v;
	
	Matrix M1(3,2);
	M1(1,1)=5.0;
	
	Matrix M2(3,2);
	M2(1,1)=3.0;
	
	Matrix M3= M1+M2;
	
	cout << "M1\n" << M1 << "\n";
	cout << "M2\n" << M2 << "\n";
	cout << "M3\n" << M3 << "\n";
	
	
	return 0;
}