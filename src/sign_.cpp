// $Source$
//--------------------------------------------------------------------------------
// sign_
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file sign_.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\sign_.hpp"

double sign_(double a,double  b){
	if (b>=0.0)
		return fabs(a);
	else
		return - fabs(a);
	

}