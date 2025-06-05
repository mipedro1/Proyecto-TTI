// $Header$
//--------------------------------------------------------------------------------
// gmst
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file gmst.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Gmst_
#define _Gmst_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Frac.hpp"
#include "../include/global.hpp"
#include <cmath>

using namespace std;

//-----------------------------------------------------------------------------------------------
// gmst(double Mjd_UT1)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Greenwich Mean Sidereal Time
 *
 *	@param [in] Mjd_UT1    Modified Julian Date UT1
 *
 *	@return double gmstime	   GMST in [rad]
 *
 */
//-----------------------------------------------------------------------------------------------
double gmst(double Mjd_UT1);

#endif