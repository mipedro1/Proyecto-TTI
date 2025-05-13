#include "..\include\MeasUpdate.hpp"

tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix& x, double z,double g,double s,Matrix& G,Matrix& P, int n){
	int m = 1;
	double Inv_W;
	
	
	Inv_W = s*s;    // Inverse weight (measurement covariance)
	

	// Kalman gain
	Matrix& ma=G*P*G.transpose()+Inv_W;
	Matrix& K = P*G.transpose()*ma.inv();
	
	// State update
	x = x + K*(z-g);
	
	// Covariance update
	P = (eye(n)-K*G)*P;
	return tie(K, x, P);
}
