// $Source$
//--------------------------------------------------------------------------------
// NutMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file NutMatrix.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\NutMatrix.hpp"

Matrix& NutMatrix (double Mjd_TT){
	// Mean obliquity of the ecliptic
	double eps = MeanObliquity (Mjd_TT);

	// Nutation in longitude and obliquity
	auto [dpsi, deps] = NutAngles (Mjd_TT);

	// Transformation from mean to true equator and equinox
	Matrix &NutMat = R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);
	return NutMat;
}

