// $Source$
//--------------------------------------------------------------------------------
// MeanObliquity
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file MeanObliquity.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\MeanObliquity.hpp"

double MeanObliquity (double Mjd_TT){
	

	double T = (Mjd_TT-MJD_J2000)/36525.0;

	return Rad *( 84381.448/3600.0-(46.8150+(0.00059-0.001813*T)*T)*T/3600.0 );

}