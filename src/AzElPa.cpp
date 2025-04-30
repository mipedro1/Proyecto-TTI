#include "..\include\AzElPa.hpp"

tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix& s){
	
	double pi2 = 2.0*M_PI;

	double rho = sqrt(s(1)*s(1)+s(2)*s(2));

	// Angles
	double Az = atan2(s(1),s(2));

	if (Az<0.0) 
		Az = Az+pi2;
	

	double El = atan ( s(3) / rho );

	// Partials
	Matrix &dAds = zeros(3);
	dAds(1,1) =  s(2)/(rho*rho);
	dAds(1,2)=-s(1)/(rho*rho);
	dAds(1,3)=0.0;
	Matrix &dEds = zeros(3);
	dEds(1,1)=-s(1)*s(3)/rho;
	dEds(1,2)=-s(2)*s(3)/rho;
	dEds(1,3)=rho;
	dEds = dEds / s.dot(s);
	
	return tie(Az, El, dAds, dEds);
	
}