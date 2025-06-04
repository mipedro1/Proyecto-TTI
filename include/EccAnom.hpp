// $Header$
//--------------------------------------------------------------------------------
// EccAnom
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file EccAnom.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _EccAnom_
#define _EccAnom_

#include "..\include\matrix.hpp"
#include <cmath>
#include <cfloat>
#include "..\include\SAT_Const.hpp"


//-----------------------------------------------------------------------------------------------
// EccAnom (double M, double e)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the eccentric anomaly for elliptic orbits
 *
 *	@param [in] M         Mean anomaly in [rad]
 *	@param [in] e         Eccentricity of the orbit [0,1]    
 *
 *	@return double E   	  Eccentric anomaly in [rad]
 *
 */
//-----------------------------------------------------------------------------------------------
double EccAnom (double M, double e);

#endif