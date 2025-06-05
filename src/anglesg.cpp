// $Source$
//--------------------------------------------------------------------------------
// anglesg
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file anglesg.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\anglesg.hpp"

tuple<Matrix&,Matrix&> anglesg ( double az1,double az2,double az3,double el1,double el2,double el3,double Mjd1,double Mjd2,double Mjd3,Matrix& Rs1,Matrix& Rs2,Matrix& Rs3 ){
	
	Matrix& L1=zeros(3,1);
	L1(1,1)=cos(el1)*sin(az1); L1(2,1)=cos(el1)*cos(az1); L1(3,1)=sin(el1);
	
	Matrix& L2=zeros(3,1);
	L2(1,1)=cos(el2)*sin(az2); L2(2,1)=cos(el2)*cos(az2); L2(3,1)=sin(el2);
	
	Matrix& L3=zeros(3,1);
	L3(1,1)=cos(el3)*sin(az3); L3(2,1)=cos(el3)*cos(az3); L3(3,1)=sin(el3);

	auto [lon1, lat1, h1] = Geodetic(Rs1);
	auto [lon2, lat2, h2] = Geodetic(Rs2);
	auto [lon3, lat3, h3] = Geodetic(Rs3);

	Matrix& M1 = LTC(lon1, lat1);
	Matrix& M2 = LTC(lon2, lat2);
	Matrix& M3 = LTC(lon3, lat3);

	// body-fixed system
	Matrix& Lb1 = M1.transpose()*L1;
	Matrix& Lb2 = M1.transpose()*L2;
	Matrix& Lb3 = M1.transpose()*L3;

	// mean of date system (J2000)
	double Mjd_UTC = Mjd1;
	auto [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(Mjd_UTC,'l');
	auto [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
	double Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	double Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

	Matrix& P = PrecMatrix(MJD_J2000,Mjd_TT);
	Matrix& N = NutMatrix(Mjd_TT);
	Matrix& T = N * P;
	Matrix& E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

	Matrix& Lm1 = E.transpose()*Lb1;
	Rs1 = E.transpose()*Rs1;

	Mjd_UTC = Mjd2;
	tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(Mjd_UTC,'l');
	tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

	P = PrecMatrix(MJD_J2000,Mjd_TT);
	N = NutMatrix(Mjd_TT);
	T = N * P;
	E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

	Matrix& Lm2 = E.transpose()*Lb2;
	Rs2 = E.transpose()*Rs2;

	Mjd_UTC = Mjd3;
	tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(Mjd_UTC,'l');
	tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;

	P = PrecMatrix(MJD_J2000,Mjd_TT);
	N = NutMatrix(Mjd_TT);
	T = N * P;
	E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

	Matrix& Lm3 = E.transpose()*Lb3;
	Rs3 = E.transpose()*Rs3;

	// geocentric inertial position
	double tau1 = (Mjd1-Mjd2)*86400.0;
	double tau3 = (Mjd3-Mjd2)*86400.0;

	double a1 = tau3/(tau3-tau1);
	double a3 =-tau1/(tau3-tau1);

	double b1 = tau3/(6.0*(tau3-tau1))*(pow((tau3-tau1),2)-pow(tau3,2));
	double b3 =-tau1/(6.0*(tau3-tau1))*(pow((tau3-tau1),2)-pow(tau1,2));
	
	Matrix& aux=zeros(3,3);
	aux.assign_column(1,Lm1); aux.assign_column(2,Lm2); aux.assign_column(3,Lm3);
	Matrix& aux2=zeros(3,3);
	aux2.assign_column(1,Rs1); aux2.assign_column(2,Rs2); aux2.assign_column(3,Rs2);
	Matrix& D = aux.inv()*aux2;

	double d1s = D(2,1)*a1-D(2,2)+D(2,3)*a3;
	double d2s = D(2,1)*b1+D(2,3)*b3;

	double Ccye = 2.0*Lm2.transpose().dot(Rs2.transpose());
	
	Matrix poly = zeros(1, 9);
	poly(1, 1) = 1.0;  // R2^8... polynomial
	poly(1, 2) = 0.0;
	poly(1, 3) = -(pow(d1s, 2) + d1s * Ccye + pow((norm(Rs2)), 2));
	poly(1, 4) = 0.0;
	poly(1, 5) = 0.0;
	poly(1, 6) = -GM_Earth * (d2s * Ccye + 2 * d1s * d2s);
	poly(1, 7) = 0.0;
	poly(1, 8) = 0.0;
	poly(1, 9) = pow(-GM_Earth, 2) * pow(d2s, 2);
	
	Matrix coef=zeros(poly.n_column);  
	Matrix zeror=zeros(poly.n_column);  
	Matrix zeroi=zeros(poly.n_column); 
		
	for (int i = 1; i <= poly.n_column; ++i) {  
		coef(1, i) = poly(1, i); 
	}
	int num_roots = real_poly_roots(&coef(1, 1), poly.n_column - 1, &zeror(1, 1), &zeroi(1, 1));

	Matrix& rootarr = zeros(num_roots, 1);
	int real_root_count = 0; 
	for (int i = 1; i <= num_roots; ++i) {
		if (fabs(zeroi(i)) < 1e-8) {  
			real_root_count++; 
			rootarr(real_root_count, 1) = zeror(i);  
		}
	}
	
	double bigr2= -99999990.0;

	for (int j=1;j<=8;j++){
		if ( rootarr(j) > bigr2  )
			bigr2= rootarr(j);
		  
	}

	double u = GM_Earth/(pow(bigr2,3));

	double C1 = a1+b1*u;
	double C2 = -1.0;
	double C3 = a3+b3*u;
	
	Matrix& cs=zeros(3,1);
	cs(1,1)=C1; cs(2,1)=C2; cs(3,1)=C3;
	Matrix& temp = D*(-1)*cs;
	double rho1 = temp(1)/(a1+b1*u);
	double rho2 = -temp(2);
	double rho3 = temp(3)/(a3+b3*u);

	double rhoold1 = rho1;
	double rhoold2 = rho2;
	double rhoold3 = rho3;

	rho2 = 99999999.9;
	int ll   = 0;
	
	Matrix& r1=zeros(3,1);
	Matrix& r2=zeros(3,1);
	Matrix& r3=zeros(3,1);
	double magr1,magr2,magr3,rdot,udot,tausqr,f1,g1,f3,g3,H1,H2,H3,G1,G2,G3,D1,D2,D3;
	Matrix& element=zeros(3,2);
	Matrix& ds=zeros(1,3);
	Matrix& cc=zeros(3,1);
	Matrix& v2 = zeros(3, 1);
	while ((fabs(rhoold2-rho2) > 1e-12) && (ll <= 0 )){
		ll = ll + 1;
		rho2 = rhoold2;
		
		r1 = Rs1+Lm1*rho1;
		r2 = Rs2+Lm2*rho2;
		r3 = Rs3+Lm3*rho3;
		
		magr1 = norm(r1);
		magr2 = norm(r2);
		magr3 = norm(r3);
		
		auto [v2, theta,theta1,copa,error] = gibbs(r1,r2,r3);
		
		if ( (error != "          ok") && (copa < M_PI/180.0) )        
			tie(v2,theta,theta1,copa,error) = hgibbs(r1,r2,r3,Mjd1,Mjd2,Mjd3);
		
		element.assign_column(1,r2);
		element.assign_column(2,v2);
		auto [p, a, e, i, Omega, omega, M] = elements (element);
		
		if ( ll <= 8 ){
			u = GM_Earth/(pow(magr2,3));
			rdot= r2.transpose().dot(v2.transpose())/magr2;
			udot= (-3.0*GM_Earth*rdot)/(pow(magr2,4));
			
			tausqr= tau1*tau1;
			f1=  1.0 - 0.5*u*tausqr -(1.0/6.0)*udot*tausqr*tau1 
				- (1.0/24.0) * u*u*tausqr*tausqr 
				- (1.0/30.0)*u*udot*tausqr*tausqr*tau1;
			g1= tau1 - (1.0/6.0)*u*tau1*tausqr - (1.0/12.0) * udot*tausqr*tausqr 
				- (1.0/120.0)*u*u*tausqr*tausqr*tau1 
				- (1.0/120.0)*u*udot*tausqr*tausqr*tausqr;
			tausqr= tau3*tau3;
			f3=  1.0 - 0.5*u*tausqr -(1.0/6.0)*udot*tausqr*tau3 
				- (1.0/24.0) * u*u*tausqr*tausqr 
				- (1.0/30.0)*u*udot*tausqr*tausqr*tau3;
			g3= tau3 - (1.0/6.0)*u*tau3*tausqr - (1.0/12.0) * udot*tausqr*tausqr 
				- (1.0/120.0)*u*u*tausqr*tausqr*tau3 
				- (1.0/120.0)*u*udot*tausqr*tausqr*tausqr;
		}else{
			
			theta  = angl( r1,r2 );
			theta1 = angl( r2,r3 );
			
			f1= 1.0 - ( (magr1*(1.0 - cos(theta)) / p ) );
			g1= ( magr1*magr2*sin(-theta) ) / sqrt( p );
			f3= 1.0 - ( (magr3*(1 - cos(theta1)) / p ) );
			g3= ( magr3*magr2*sin(theta1) ) / sqrt( p );
		}
		
		C1 = g3/(f1*g3-f3*g1);
		C2 = -1.0;
		C3 =-g1/(f1*g3-f3*g1);
		
		H1 = GM_Earth*tau3/12.0;
		H3 =-GM_Earth*tau1/12.0;
		H2 = H1-H3;
		
		G1 = -tau3/(tau1*(tau3-tau1));
		G3 = -tau1/(tau3*(tau3-tau1));
		G2 = G1-G3;
		
		D1 = G1+H1/(pow(magr1,3));
		D2 = G2+H2/(pow(magr2,3));
		D3 = G3+H3/(pow(magr3,3));
		ds(1,1)=D1; ds(1,2)=D2;ds(1,3)=D3;
		cc(1,1)=C1; cc(2,1)=C2;cc(3,1)=C3;
		temp = ds*(-1)*cc;
		rhoold1 = temp(1)/(a1+b1*u);
		rhoold2 = -temp(1);
		rhoold3 = temp(1)/(a3+b3*u);
		
		r1 = Rs1+Lm1*rhoold1;
		r2 = Rs2+Lm2*rhoold2;
		r3 = Rs3+Lm3*rhoold3;
		
	}

	r1 = Rs1+Lm1*rho1;
	r2 = Rs2+Lm2*rho2;
	r3 = Rs3+Lm3*rho3;
	
	return tie(r2,v2);
}