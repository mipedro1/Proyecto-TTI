// $Header$
//--------------------------------------------------------------------------------
// NutAngles
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file NutAngles.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _NutAngles_
#define _NutAngles_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>

using namespace std;


//-----------------------------------------------------------------------------------------------
// NutAngles (double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Nutation in longitude and obliquity
 *
 *	@param [in] Mjd_TT     Modified Julian Date (Terrestrial Time)
 *
 *	@return double dpsi , Nutation Angle
 *	@return double deps , Nutation Angle
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<double,double> NutAngles (double Mjd_TT);

#endif