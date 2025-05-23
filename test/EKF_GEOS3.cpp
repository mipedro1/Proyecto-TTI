#include "../include/matrix.hpp"
#include "../include/R_x.hpp"
#include "../include/global.hpp"
#include "..\include\Position.hpp"
#include "..\include\DEInteg.hpp"
#include "..\include\LTC.hpp"
#include "..\include\timediff.hpp"
#include "..\include\IERS.hpp"
#include "..\include\gmst.hpp"
#include "..\include\R_x.hpp"
#include "..\include\TimeUpdate.hpp"
#include "..\include\AzElPa.hpp"
#include "..\include\MeasUpdate.hpp"
#include "..\include\Accel.hpp"
#include "..\include\VarEqn.hpp"
#include <cstdio>
#include <cmath>
#include <tuple>

using namespace std;


int main(){
	eop19620101(21413);
	
	GGM03S(181);
	
	DE430Coeff(2285,1020);
	
	GEOS3(46);
	
	
	double sigma_range,sigma_az,sigma_el,lat,lon,alt,Mjd1,Mjd2,Mjd3,theta,Dist,Azim,Elev,x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC,UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
	
	sigma_range = 92.5;          // [m]
	sigma_az = 0.0224*Rad; // [rad]
	sigma_el = 0.0139*Rad; // [rad]

	// Kaena Point station
	lat = Rad*21.5748;     // [rad]
	lon = Rad*(-158.2706); // [rad]
	alt = 300.20;                // [m]

	Matrix& Rs = Position(lon, lat, alt).transpose();

	Mjd1 = obs(1,1);
	Mjd2 = obs(9,1);
	Mjd3 = obs(18,1);

	
	Matrix &r2=zeros(3),&v2=zeros(3);
	r2(1)=6221397.62857869; r2(2)=2867713.77965738; r2(3)=3006155.98509949;
	v2(1)=4645.04725161806; v2(2)=-2752.21591588204; v2(3)=-7507.99940987031;
	//auto [r2,v2] = anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);
	// [r2,v2] = anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
	//                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs);

	Matrix& Y0_apr = zeros(6,1);
	Y0_apr(1,1)=r2(1); Y0_apr(2,1)=r2(2); Y0_apr(3,1)=r2(3);
	Y0_apr(4,1)=v2(1); Y0_apr(5,1)=v2(2); Y0_apr(6,1)=v2(3);

	double Mjd0 = Mjday(1995,1,29,02,38,0);

	double Mjd_UTC = obs(9,1);
	double Mjd_TT,Mjd_UT1;
	AuxParam.Mjd_UTC = Mjd_UTC;
	AuxParam.n      = 20;
	AuxParam.m      = 20;
	AuxParam.sun     = 1;
	AuxParam.moon    = 1;
	AuxParam.planets = 1;

	int n_eqn  = 6;
	Matrix& Y = DEInteg(Accel,0.0,-(obs(9,1)-Mjd0)*86400.0,1e-13,1e-6,6,Y0_apr);
	Matrix& P = zeros(6,6);
	  
	for (int i=1;i<=3;i++){
		P(i,i)=1e8;
	}
	for (int i=4;i<=6;i++){
		P(i,i)=1e3;
	}
	Matrix& LT = LTC(lon,lat);
	Matrix& yPhi = zeros(42,1);
	Matrix& Phi  = zeros(6,6);
	
	// Measurement loop
	int t = 0;
	int nobs = 46;
	Matrix& Y_old=zeros(6,1);
	Matrix& U=zeros(6,1);
	Matrix& r=zeros(6,1);
	Matrix& s=zeros(6,1);
	Matrix& dAdY=zeros(6,1);
	Matrix& dEdY=zeros(6,1);
	Matrix& dDds=zeros(6,1);
	Matrix& dDdY=zeros(6,1);
	Matrix& K=zeros(6,1);
	Matrix& dAds=zeros(6,1);
	Matrix& dEds=zeros(6,1);
	int t_old;
	for (int i=1;i<=nobs;i++){
		// Previous step
		t_old = t;
		Y_old = Y;
		
		// Time increment and propagation
		Mjd_UTC = obs(i,1);                       // Modified Julian Date
		t       = (Mjd_UTC-Mjd0)*86400.0;         // Time since epoch [s]
		
		tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(Mjd_UTC,'l');
		tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
		Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
		Mjd_UT1 = Mjd_TT + (UT1_UTC-TT_UTC)/86400.0;
		AuxParam.Mjd_UTC = Mjd_UTC;
		AuxParam.Mjd_TT = Mjd_TT;
			
		for (int ii=1;ii<=6;ii++){
			yPhi(ii) = Y_old(ii);
			for (int j=1;j<=6;j++){
				if (ii==j) 
					yPhi(6*j+ii) = 1; 
				else
					yPhi(6*j+ii) = 0;
				
			}
		}
		
		yPhi = DEInteg (VarEqn,0,t-t_old,1e-13,1e-6,42,yPhi);
		
		// Extract state transition matrices
		for (int j=1;j<=6;j++){
			Phi.assign_column(j, yPhi.transpose().extract_vector(6*j+1, 6*j+6).transpose());
		}
		
		Y = DEInteg (Accel,0,t-t_old,1e-13,1e-6,6,Y_old);
		
		// Topocentric coordinates
		theta = gmst(Mjd_UT1);                    // Earth rotation
		U = R_z(theta);
		r = Y.transpose().extract_vector(1,3);
		s = LT*(U*r.transpose()-Rs);                          // Topocentric position [m]
		
		// Time update
		P = TimeUpdate(P, Phi);
			
		// Azimuth and partials
		tie(Azim, Elev, dAds, dEds) = AzElPa(s.transpose());     // Azimuth, Elevation
		dAdY = (dAds * LT * U).union_vector(zeros(1, 3));
		
		
		// Measurement update
		tie (K, Y, P) = MeasUpdate ( Y, obs(i,2), Azim, sigma_az, dAdY, P, 6 );
		// Elevation and partials
		r = Y.transpose().extract_vector(1, 3);
		s = LT*(U*r.transpose()-Rs);                          // Topocentric position [m]
		tie(Azim, Elev, dAds, dEds) = AzElPa(s.transpose());     // Azimuth, Elevation
		dEdY =(dEds*LT*U).union_vector(zeros(1,3));
		// Measurement update
		tie(K, Y, P) = MeasUpdate ( Y, obs(i,3), Elev, sigma_el, dEdY, P, 6 );
		
		// Range and partials
		r = Y.transpose().extract_vector(1, 3);
		s = LT*(U*r.transpose()-Rs);                          // Topocentric position [m]
		Dist = norm(s); dDds = (s/Dist).transpose();         // Range
		dDdY = (dDds*LT*U).union_vector(zeros(1,3));
		
		// Measurement update
		tie(K, Y, P) = MeasUpdate ( Y, obs(i,4), Dist, sigma_range, dDdY, P, 6 );
	}
	
	tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC) = IERS(obs(46,1),'l');
	tie(UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC) = timediff(UT1_UTC,TAI_UTC);
	Mjd_TT = Mjd_UTC + TT_UTC/86400.0;
	AuxParam.Mjd_UTC = Mjd_UTC;
	AuxParam.Mjd_TT = Mjd_TT;

	Matrix& Y0 = DEInteg (Accel,0,-(obs(46,1)-obs(1,1))*86400.0,1e-13,1e-6,6,Y);
	
	Matrix& Y_true=zeros(6,1);
	Y_true(1,1)=5753.173e3;
	Y_true(2,1)=2673.361e3;
	Y_true(3,1)=3440.304e3;
	Y_true(4,1)=4.324207e3;
	Y_true(5,1)=-1.924299e3;
	Y_true(6,1)=-5.728216e3;

	printf("\nError of Position Estimation\n");
	printf("dX%10.1lf [m]\n", Y0(1) - Y_true(1));
	printf("dY%10.1lf [m]\n", Y0(2) - Y_true(2));
	printf("dZ%10.1lf [m]\n", Y0(3) - Y_true(3));
	printf("\nError of Velocity Estimation\n");
	printf("dVx%8.1lf [m/s]\n", Y0(4) - Y_true(4));
	printf("dVy%8.1lf [m/s]\n", Y0(5) - Y_true(5));
	printf("dVz%8.1lf [m/s]\n", Y0(6) - Y_true(6));
	
	
	return 0;
}