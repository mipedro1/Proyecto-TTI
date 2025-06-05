// $Header$
//--------------------------------------------------------------------------------
// AccelPointMass
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file AccelPointMass.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _AccelPointMass_
#define _AccelPointMass_

#include "..\include\matrix.hpp"
#include <cmath>


//-----------------------------------------------------------------------------------------------
// AccelPointMass(Matrix& r, Matrix& s,double GM)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the perturbational acceleration due to a point
 *		   mass
 *
 *	@param [in] r           Satellite position vector 
 *	@param [in] s           Point mass position vector
 *	@param [in] GM          Gravitational coefficient of point mass
 *
 *	@return Matrix& a           Acceleration (a=d^2r/dt^2)
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& AccelPointMass(Matrix& r, Matrix& s,double GM);

#endif