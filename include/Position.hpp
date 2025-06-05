// $Header$
//--------------------------------------------------------------------------------
// Position
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file Position.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Position_
#define _Position_

#include "..\include\matrix.hpp"
#include <cmath>
#include "..\include\SAT_Const.hpp"


//-----------------------------------------------------------------------------------------------
// Position(double lon,double  lat, double h)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Position vector (r [m]) from geodetic coordinates (Longitude [rad],
 *		   latitude [rad], altitude [m])
 *
 *	@param [in]  lon
 *	@param [in]  lat   
 *  @param [in]  h 			 
 *
 *	@return Matrix& r
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& Position(double lon,double  lat, double h);

#endif