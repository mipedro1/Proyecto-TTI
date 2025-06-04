// $Header$
//--------------------------------------------------------------------------------
// LTC
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file LTC.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _LTC_
#define _LTC_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include <cmath>

using namespace std;


//-----------------------------------------------------------------------------------------------
// LTC(double lon,double lat)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Transformation from Greenwich meridian system to 
 *	       local tangent coordinates
 *
 *	@param [in] lon      -Geodetic East longitude [rad]
 *	@param [in] lat      -Geodetic latitude [rad]
 *
 *	@return Matrix& M        -Rotation matrix from the Earth equator and Greenwich meridian
 *							 to the local tangent (East-North-Zenith) coordinate system
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& LTC(double lon,double lat);

#endif