// $Source$
//--------------------------------------------------------------------------------
// EqnEquinox
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file EqnEquinox.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\EqnEquinox.hpp"

double EqnEquinox (double Mjd_TT){
	// Nutation in longitude and obliquity
	auto [dpsi, deps] = NutAngles (Mjd_TT);

	// Equation of the equinoxes
	return dpsi * cos ( MeanObliquity(Mjd_TT) );

}