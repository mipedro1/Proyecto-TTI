// $Header$
//--------------------------------------------------------------------------------
// gast
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file gast.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Gast_
#define _Gast_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include "..\include\EqnEquinox.hpp"
#include "..\include\gmst.hpp"
#include <cmath>

using namespace std;


//-----------------------------------------------------------------------------------------------
// gast (double Mjd_UT1)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Greenwich Apparent Sidereal Time
 *
 *	@param [in] Mjd_UT1   Modified Julian Date UT1
 *
 *	@return double gstime    GAST in [rad]
 *
 */
//-----------------------------------------------------------------------------------------------
double gast (double Mjd_UT1);

#endif