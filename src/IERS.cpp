#include "..\include\IERS.hpp"

tuple<double,double,double,double,double,double,double,double,double> IERS(double Mjd_UTC,char interp='n'){
	
	double mjd, mfme, fixf, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
	int i;
	if (interp =='l'){
		// linear interpolation
		mjd = (floor(Mjd_UTC));
		i = -1;  
		for (int j = 1; j <= eopdata.n_column; j++) {
			if (mjd == eopdata(4, j)) {
				i = j;  
				break;   
			}
		}
		Matrix preeopdata = eopdata.extract_column(i);      
		Matrix nexteopdata = eopdata.extract_column(i+1);   
		 mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
		 fixf = mfme/1440;
		// Setting of IERS Earth rotation parameters
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		 x_pole  = preeopdata(5)+(nexteopdata(5)-preeopdata(5))*fixf;
		 y_pole  = preeopdata(6)+(nexteopdata(6)-preeopdata(6))*fixf;
		 UT1_UTC = preeopdata(7)+(nexteopdata(7)-preeopdata(7))*fixf;
		 LOD     = preeopdata(8)+(nexteopdata(8)-preeopdata(8))*fixf;
		 dpsi    = preeopdata(9)+(nexteopdata(9)-preeopdata(9))*fixf;
		 deps    = preeopdata(10)+(nexteopdata(10)-preeopdata(10))*fixf;
		 dx_pole = preeopdata(11)+(nexteopdata(11)-preeopdata(11))*fixf;
		 dy_pole = preeopdata(12)+(nexteopdata(12)-preeopdata(12))*fixf;
		 TAI_UTC = preeopdata(13);
		
		x_pole  = x_pole/Arcs;  // Pole coordinate [rad]
		y_pole  = y_pole/Arcs;  // Pole coordinate [rad]
		dpsi    = dpsi/Arcs;
		deps    = deps/Arcs;
		dx_pole = dx_pole/Arcs; // Pole coordinate [rad]
		dy_pole = dy_pole/Arcs; // Pole coordinate [rad]
	}else if (interp =='n') {   
		mjd = (floor(Mjd_UTC));
		int i = -1;
		for (int j = 0; j < eopdata.n_column; j++) {
			if (eopdata(4, j) == mjd) {
				i = j;
				break;  
			}
		}
		Matrix eopdata_col = eopdata.extract_column(i);
		// Setting of IERS Earth rotation parameters
		// (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
		x_pole  = eopdata(5)/Arcs;  // Pole coordinate [rad]
		y_pole  = eopdata(6)/Arcs;  // Pole coordinate [rad]
		UT1_UTC = eopdata(7);             // UT1-UTC time difference [s]
		LOD     = eopdata(8);             // Length of day [s]
		dpsi    = eopdata(9)/Arcs;
		deps    = eopdata(10)/Arcs;
		dx_pole = eopdata(11)/Arcs; // Pole coordinate [rad]
		dy_pole = eopdata(12)/Arcs; // Pole coordinate [rad]
		TAI_UTC = eopdata(13);            // TAI-UTC time difference [s]
	}
	
	return tie(x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC);
}