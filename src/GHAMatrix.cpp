// $Source$
//--------------------------------------------------------------------------------
// GHAMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file GHAMatrix.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\GHAMatrix.hpp"

Matrix& GHAMatrix (double Mjd_UT1){
	Matrix& GHAmat = R_z( gast(Mjd_UT1) );
	return GHAmat;
}