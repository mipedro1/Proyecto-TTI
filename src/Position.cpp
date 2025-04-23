#include "..\include\Position.hpp"


Matrix Position(double lon,double  lat, double h){
	
	double R_equ = R_Earth;
	double f     = f_Earth;

	double e2     = f*(2.0-f);   // Square of eccentricity
	double CosLat = cos(lat);    // (Co)sine of geodetic latitude
	double SinLat = sin(lat);

	// Position vector 
	double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);
	
	Matrix r(3);
	r(1) =  (         N+h)*CosLat*cos(lon);
	r(2) =  (         N+h)*CosLat*sin(lon);
	r(3) =  ((1.0-e2)*N+h)*SinLat;
	
	return r;
}