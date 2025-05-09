#include "..\include\JPL_Eph_DE430.hpp"

tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&> JPL_Eph_DE430(double Mjd_TDB){
	double JD,t1,dt,Mjd0;
	JD = Mjd_TDB + 2400000.5;

	int i;
	for(i=1;i<=PC.n_row;i++){
		if(PC(i,1)<=JD && JD<=PC(i,2)){
			break;
		}
	}
	Matrix PCtemp = PC.extract_row(i);

	t1 = PCtemp(1)-2400000.5; // MJD at start of interval

	dt = Mjd_TDB - t1;
	
	Matrix& temp=zeros(4);
	int j;
	for(j=0;j<4;j++){
		temp(j+1)=231+13*j;
	}
	Matrix Cx_Earth,Cy_Earth,Cz_Earth,Cx,Cy,Cz;
	Cx_Earth = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Earth = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Earth = PCtemp.extract_vector(temp(3), temp(4) - 1);
	temp = temp+39;
	Cx = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz = PCtemp.extract_vector(temp(3), temp(4) - 1);
	Cx_Earth = Cx_Earth.union_vector(Cx);
	Cy_Earth = Cy_Earth.union_vector(Cy);
	Cz_Earth = Cz_Earth.union_vector(Cz);    
	if (0<=dt && dt<=16){
		j=0;
		Mjd0 = t1;
	}else if(16<dt && dt<=32){
		j=1;
		Mjd0 = t1+16*j;
	}
	Matrix& r_Earth = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+16, Cx_Earth.extract_vector(13 * j + 1, 13 * j + 13),
						 Cy_Earth.extract_vector(13 * j + 1, 13 * j + 13), Cz_Earth.extract_vector(13 * j + 1, 13 * j + 13)).transpose()*1e3;
	

	for (int j = 0; j < 4; j++) {
		temp(j+1) = 441 + 13 * j;  
	}
	Matrix Cx_Moon,Cy_Moon,Cz_Moon;
	Cx_Moon = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Moon = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Moon = PCtemp.extract_vector(temp(3), temp(4) - 1);
	for (int i = 1; i <= 7; i++) {
		temp = temp+39;
		Cx = PCtemp.extract_vector(temp(1), temp(2) - 1);
		Cy = PCtemp.extract_vector(temp(2), temp(3) - 1);
		Cz = PCtemp.extract_vector(temp(3), temp(4) - 1);
		Cx_Moon = Cx_Moon.union_vector(Cx);
		Cy_Moon = Cy_Moon.union_vector(Cy);
		Cz_Moon = Cz_Moon.union_vector(Cz); 
	}
	if (0<=dt && dt<=4){
		j=0;
		Mjd0 = t1;
	}else if(4<dt && dt<=8){
		j=1;
		Mjd0 = t1+4*j;
	}else if(8<dt && dt<=12){
		j=2;
		Mjd0 = t1+4*j;
	}else if(12<dt && dt<=16){
		j=3;
		Mjd0 = t1+4*j;
	}else if(16<dt && dt<=20){
		j=4;
		Mjd0 = t1+4*j;
	}else if(20<dt && dt<=24){
		j=5;
		Mjd0 = t1+4*j;
	}else if(24<dt && dt<=28){
		j=6;
		Mjd0 = t1+4*j;
	}else if(28<dt && dt<=32){
		j=7;
		Mjd0 = t1+4*j;
	}
	Matrix& r_Moon = Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, Cx_Moon.extract_vector(13 * j + 1, 13 * j + 13),
						Cy_Moon.extract_vector(13 * j + 1, 13 * j + 13), Cz_Moon.extract_vector(13 * j + 1, 13 * j + 13)).transpose()*1e3;

	
	for (int j = 0; j < 4; j++) {
		temp(j + 1) = 753 + 11 * j;
	}
	Matrix Cx_Sun, Cy_Sun, Cz_Sun;
	Cx_Sun = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Sun = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Sun = PCtemp.extract_vector(temp(3), temp(4) - 1);
	temp = temp + 33;
	Cx = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz = PCtemp.extract_vector(temp(3), temp(4) - 1);  
	Cx_Sun = Cx_Sun.union_vector(Cx);
	Cy_Sun = Cy_Sun.union_vector(Cy);
	Cz_Sun = Cz_Sun.union_vector(Cz);
	if (0<=dt && dt<=16){
		j=0;
		Mjd0 = t1;
	}else if(16<dt && dt<=32){
		j=1;
		Mjd0 = t1+16*j;
	}
	Matrix& r_Sun = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, Cx_Sun.extract_vector(11 * j + 1, 11 * j + 11),
					   Cy_Sun.extract_vector(11 * j + 1, 11 * j + 11), Cz_Sun.extract_vector(11 * j + 1, 11 * j + 11)).transpose()*1e3;

	
	for (int j = 0; j < 4; j++) {
		temp(j + 1) = 3 + 14 * j;
	}
	Matrix Cx_Mercury, Cy_Mercury, Cz_Mercury;
	Cx_Mercury = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Mercury = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Mercury = PCtemp.extract_vector(temp(3), temp(4) - 1);
	for (int i = 1; i < 3; i++) {
		temp = temp + 42;
		Cx = PCtemp.extract_vector(temp(1), temp(2) - 1);
		Cy = PCtemp.extract_vector(temp(2), temp(3) - 1);
		Cz = PCtemp.extract_vector(temp(3), temp(4) - 1);
		
		Cx_Mercury = Cx_Mercury.union_vector(Cx);
		Cy_Mercury = Cy_Mercury.union_vector(Cy);
		Cz_Mercury = Cz_Mercury.union_vector(Cz);
	}
	if (0<=dt && dt<=8){
		j=0;
		Mjd0 = t1;
	}else if(8<dt && dt<=16){
		j=1;
		Mjd0 = t1+8*j;
	}else if (16<dt && dt<=24){
		j=2;
		Mjd0 = t1+8*j;
	}else if(24<dt && dt<=32){
		j=3;
		Mjd0 = t1+8*j;
	}
	Matrix& r_Mercury = Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, Cx_Mercury.extract_vector(14 * j + 1, 14 * j + 14),
						   Cy_Mercury.extract_vector(14 * j + 1, 14 * j + 14), Cz_Mercury.extract_vector(14 * j + 1, 14 * j + 14)).transpose()*1e3;

	
	for (int j = 0; j < 4; j++) {
		temp(j + 1) = 171 + 10 * j;
	}
	Matrix Cx_Venus, Cy_Venus, Cz_Venus;
	Cx_Venus = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Venus = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Venus = PCtemp.extract_vector(temp(3), temp(4) - 1);
	temp = temp+30;
	Cx = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz = PCtemp.extract_vector(temp(3), temp(4) - 1);
	Cx_Venus = Cx_Venus.union_vector(Cx);
	Cy_Venus = Cy_Venus.union_vector(Cy);
	Cz_Venus = Cz_Venus.union_vector(Cz);
	if (0<=dt && dt<=16){
		j=0;
		Mjd0 = t1;
	}else if(16<dt && dt<=32){
		j=1;
		Mjd0 = t1+16*j;
	}
	Matrix & r_Venus = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, Cx_Venus.extract_vector(10 * j + 1, 10 * j + 10),
						 Cy_Venus.extract_vector(10 * j + 1, 10 * j + 10),Cz_Venus.extract_vector(10 * j + 1, 10 * j + 10)).transpose()*1e3;

	
	for (int j = 0; j < 4; j++) {
		temp(j + 1) = 309 + 11 * j;
	}
	Matrix Cx_Mars, Cy_Mars, Cz_Mars;
	Cx_Mars = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Mars = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Mars = PCtemp.extract_vector(temp(3), temp(4) - 1);
	j=0;
	Mjd0 = t1;
	Matrix& r_Mars = Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, Cx_Mars.extract_vector(11 * j + 1, 11 * j + 11),
						Cy_Mars.extract_vector(11 * j + 1, 11 * j + 11), Cz_Mars.extract_vector(11 * j + 1, 11 * j + 11)).transpose()*1e3;

	
	for(int j = 0; j < 4; j++) {
		temp(j + 1) = 342 + 8 * j;
	}
	Matrix Cx_Jupiter, Cy_Jupiter, Cz_Jupiter;
	Cx_Jupiter = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Jupiter = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Jupiter = PCtemp.extract_vector(temp(3), temp(4) - 1);
	j=0;
	Mjd0 = t1;
	Matrix& r_Jupiter = Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, Cx_Jupiter.extract_vector(8 * j + 1, 8 * j + 8),
						   Cy_Jupiter.extract_vector(8 * j + 1, 8 * j + 8), Cz_Jupiter.extract_vector(8 * j + 1, 8 * j + 8)).transpose()*1e3;

	
	for (int j = 0; j < 4; j++) {
		temp(j + 1) = 366 + 7 * j;
	}
	Matrix Cx_Saturn, Cy_Saturn, Cz_Saturn;
	Cx_Saturn = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Saturn = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Saturn = PCtemp.extract_vector(temp(3), temp(4) - 1);
	j=0;
	Mjd0 = t1;
	Matrix& r_Saturn = Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, Cx_Saturn.extract_vector(7*j+1, 7*j+7),
						  Cy_Saturn.extract_vector(7*j+1, 7*j+7), Cz_Saturn.extract_vector(7*j+1, 7*j+7)).transpose()*1e3;
	
	
	
	for(int j = 0; j < 4; j++) {
		temp(j + 1) = 387 + 6 * j;
	}
	Matrix Cx_Uranus, Cy_Uranus, Cz_Uranus;
	Cx_Uranus = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Uranus = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Uranus = PCtemp.extract_vector(temp(3), temp(4) - 1);
	j=0;
	Mjd0 = t1;
	Matrix& r_Uranus = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Uranus.extract_vector(6 * j + 1, 6 * j + 6),
						  Cy_Uranus.extract_vector(6 * j + 1, 6 * j + 6), Cz_Uranus.extract_vector(6 * j + 1, 6 * j + 6)).transpose()*1e3;

	
	for(int j = 0; j < 4; j++) {
		temp(j + 1) = 405 + 6 * j;
	}
	Matrix Cx_Neptune, Cy_Neptune, Cz_Neptune;
	Cx_Neptune = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Neptune = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Neptune = PCtemp.extract_vector(temp(3), temp(4) - 1);
	j=0;
	Mjd0 = t1;
	Matrix& r_Neptune = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Neptune.extract_vector(6 * j + 1, 6 * j + 6),
						   Cy_Neptune.extract_vector(6 * j + 1, 6 * j + 6), Cz_Neptune.extract_vector(6 * j + 1, 6 * j + 6)).transpose()*1e3;

	
	for (int j = 0; j < 4; j++) {
		temp(j + 1) = 423 + 6 * j;
	}
	Matrix Cx_Pluto, Cy_Pluto, Cz_Pluto;
	Cx_Pluto = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Pluto = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Pluto = PCtemp.extract_vector(temp(3), temp(4) - 1);
	j=0;
	Mjd0 = t1;
	Matrix& r_Pluto = Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, Cx_Pluto.extract_vector(6 * j + 1, 6 * j + 6),
						 Cy_Pluto.extract_vector(6 * j + 1, 6 * j + 6), Cz_Pluto.extract_vector(6 * j + 1, 6 * j + 6)).transpose()*1e3;

	
	for (int j = 0; j < 4; j++) {
		temp(j + 1) = 819 + 10 * j;
	}

	Matrix Cx_Nutations, Cy_Nutations;
	Cx_Nutations = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Nutations = PCtemp.extract_vector(temp(2), temp(3) - 1);
	for (int i = 0; i < 3; i++) {
		temp = temp+20;
		Cx = PCtemp.extract_vector(temp(1), temp(2) - 1);
		Cy = PCtemp.extract_vector(temp(2), temp(3) - 1);
		Cx_Nutations = Cx_Nutations.union_vector(Cx);
		Cy_Nutations = Cy_Nutations.union_vector(Cy);
	}
	if (0<=dt && dt<=8){
		j=0;
		Mjd0 = t1;
	}else if(8<dt && dt<=16){
		j=1;
		Mjd0 = t1+8*j;
	}else if (16<dt && dt<=24){
		j=2;
		Mjd0 = t1+8*j;
	}else if(24<dt && dt<=32){
		j=3;
		Mjd0 = t1+8*j;
	}
	Matrix& Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Nutations.extract_vector(10*j + 1, 10*j + 10),
					   Cy_Nutations.extract_vector(10*j + 1, 10*j + 10),zeros(10,1)).transpose();
    

	for (int j = 0; j < 4; j++) {
		temp(j + 1) = 899 + 10 * j;
	}
	Matrix Cx_Librations,Cy_Librations,Cz_Librations;
	Cx_Librations = PCtemp.extract_vector(temp(1), temp(2) - 1);
	Cy_Librations = PCtemp.extract_vector(temp(2), temp(3) - 1);
	Cz_Librations = PCtemp.extract_vector(temp(3), temp(4) - 1);
	for (int i = 0; i < 3; i++) {
		temp = temp+30;
		Cx = PCtemp.extract_vector(temp(1), temp(2) - 1);
		Cy = PCtemp.extract_vector(temp(2), temp(3) - 1);
		Cz = PCtemp.extract_vector(temp(3), temp(4) - 1);
		Cx_Librations = Cx_Librations.union_vector(Cx);
		Cy_Librations = Cy_Librations.union_vector(Cy);
		Cz_Librations = Cz_Librations.union_vector(Cz);    
	}
	if (0<=dt && dt<=8){
		j=0;
		Mjd0 = t1;
	}else if(8<dt && dt<=16){
		j=1;
		Mjd0 = t1+8*j;
	}else if (16<dt && dt<=24){
		j=2;
		Mjd0 = t1+8*j;
	}else if(24<dt && dt<=32){
		j=3;
		Mjd0 = t1+8*j;
	}
	Matrix& Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, Cx_Librations.extract_vector(10*j+1, 10*j+10),
						Cy_Librations.extract_vector(10*j+1, 10*j+10), Cz_Librations.extract_vector(10*j+1, 10*j+10)).transpose();
	
	double EMRAT = 81.30056907419062; // DE430
	double EMRAT1 = 1/(1+EMRAT);
	r_Earth = r_Earth-r_Moon*EMRAT1;
	r_Mercury = r_Earth*(-1)+r_Mercury;
	r_Venus = r_Earth*(-1)+r_Venus;
	r_Mars = r_Earth*(-1)+r_Mars;
	r_Jupiter = r_Earth*(-1)+r_Jupiter;
	r_Saturn = r_Earth*(-1)+r_Saturn;
	r_Uranus = r_Earth*(-1)+r_Uranus;
	r_Neptune = r_Earth*(-1)+r_Neptune;
	r_Pluto = r_Earth*(-1)+r_Pluto;
	r_Sun = r_Earth*(-1)+r_Sun;

	return tie(r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus,r_Neptune,r_Pluto,r_Moon,r_Sun);
}
