// $Source$
//--------------------------------------------------------------------------------
// PoleMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file PoleMatrix.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\PoleMatrix.hpp"

Matrix& PoleMatrix (double xp,double yp){
	
	Matrix& PoleMat = R_y(-xp) * R_x(-yp);
	return PoleMat;
}