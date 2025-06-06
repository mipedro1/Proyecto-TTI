// $Source$
//--------------------------------------------------------------------------------
// AccelHarmonic
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file AccelHarmonic.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//-------------------------------------------------------------------------------- 
#include "..\include\AccelHarmonic.hpp"

 Matrix& AccelHarmonic(Matrix& r, Matrix& E, double n_max,double m_max){
	double d,latgc,lon,dUdr,dUdlatgc,dUdlon,q3,q2,q1,b1,b2,b3,r2xy,ax,ay,az;
	
	const double r_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
	const double gm    = 398600.4415e9; // [m^3/s^2]; GGM03S

	// Body-fixed position
	Matrix r_bf = E * r;
	
	// Auxiliary quantities
	d = norm(r_bf);                     // distance
	
	latgc = asin(r_bf(3)/d);
	
	lon = atan2(r_bf(2),r_bf(1));
	

	
	auto [pnm, dpnm] = Legendre(n_max,m_max,latgc);
	
	
	dUdr = 0.0;
	dUdlatgc = 0.0;
	dUdlon = 0.0;
	q3 = 0.0; q2 = q3; q1 = q2;
	int n;
	int m;
	for (n=0;n<=n_max;n++){
		b1 = (-gm/pow(d,2))*pow((r_ref/d),n)*(n+1);
		b2 =  (gm/d)*pow((r_ref/d),n);
		b3 =  (gm/d)*pow((r_ref/d),n);
		for (int m=0;m<=m_max;m++){
			q1 = q1 + pnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
			q2 = q2 + dpnm(n+1,m+1)*(Cnm(n+1,m+1)*cos(m*lon)+Snm(n+1,m+1)*sin(m*lon));
			q3 = q3 + m*pnm(n+1,m+1)*(Snm(n+1,m+1)*cos(m*lon)-Cnm(n+1,m+1)*sin(m*lon));
			
		}
		dUdr     = dUdr     + q1*b1;
		dUdlatgc = dUdlatgc + q2*b2;
		dUdlon   = dUdlon   + q3*b3;
		q3 = 0; q2 = q3; q1 = q2;
	}

	// Body-fixed acceleration
	r2xy = pow(r_bf(1),2)+pow(r_bf(2),2);

	ax = (1.0/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1.0/r2xy*dUdlon)*r_bf(2);
	ay = (1.0/d*dUdr-r_bf(3)/(pow(d,2)*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1.0/r2xy*dUdlon)*r_bf(1);
	az =  1.0/d*dUdr*r_bf(3)+sqrt(r2xy)/pow(d,2)*dUdlatgc;
	
	
	
	Matrix a_bf(3,1);
	a_bf(1)=ax; a_bf(2)=ay; a_bf(3)=az;
	// Inertial acceleration 
	Matrix& a = E.transpose()*a_bf;
	return a;
 }