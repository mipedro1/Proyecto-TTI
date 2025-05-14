#include "..\include\G_AccelHarmonic.hpp"

Matrix& G_AccelHarmonic( Matrix& r,Matrix& U, int n_max,int  m_max ){
	double d = 1.0;   // Position increment [m]

	Matrix& G = zeros(3,3);
	Matrix& dr = zeros(3,1);
	Matrix& da=zeros(3);
	Matrix& zero_value=zeros(3,1);
	// Gradient
	for (int i=1;i<=3;i++){
		// Set offset in i-th component of the position vector
		dr.assign_column(1, zero_value);
		dr(i,1) = d;
		// Acceleration difference
		da = AccelHarmonic ( r+dr/2,U, n_max, m_max ) - 
			 AccelHarmonic ( r-dr/2,U, n_max, m_max );
		// Derivative with respect to i-th axis
		G.assign_column(i, da / d); 
	}
	return G;
}