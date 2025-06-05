// $Source$
//--------------------------------------------------------------------------------
// AccelPointMass
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file AccelPointMass.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\AccelPointMass.hpp"
Matrix& AccelPointMass(Matrix& r, Matrix& s,double GM){
	
	// Relative position vector of satellite w.r.t. point mass 
	Matrix d = r - s;
	
	
	// Acceleration 
	Matrix& a = ( d/pow(norm(d),3) + s/pow(norm(s),3) ) * (-GM);
	
	
	return a;
	
	
	



}