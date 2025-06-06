// $Source$
//--------------------------------------------------------------------------------
// Position
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file Position.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\Position.hpp"


Matrix& Position(double lon,double  lat, double h){
	
	double R_equ = R_Earth;
	double f     = f_Earth;

	double e2     = f*(2.0-f);   // Square of eccentricity
	double CosLat = cos(lat);    // (Co)sine of geodetic latitude
	double SinLat = sin(lat);

	// Position vector 
	double N = R_equ / sqrt(1.0-e2*SinLat*SinLat);
	
	Matrix& r=zeros(3);
	r(1) =  (         N+h)*CosLat*cos(lon);
	r(2) =  (         N+h)*CosLat*sin(lon);
	r(3) =  ((1.0-e2)*N+h)*SinLat;
	
	return r;
}