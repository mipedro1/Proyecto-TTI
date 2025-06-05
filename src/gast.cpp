// $Source$
//--------------------------------------------------------------------------------
// gast
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file gast.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\gast.hpp"

double gast (double Mjd_UT1){
	double gstime= fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), pi2 );
	
	if(gstime<0){
		gstime+=pi2;
	}
	return gstime;
}