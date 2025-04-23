#include "../include/matrix.hpp"
#include "../include/R_x.hpp"
#include "../include/global.hpp"

int main(){
	eop19620101(4); // c=21413
	cout<< "eopdata\n"<<eopdata<<"\n";
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