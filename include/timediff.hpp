// $Header$
//--------------------------------------------------------------------------------
// timediff
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file timediff.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Timediff_
#define _Timediff_


#include <cmath>
#include <tuple>

using namespace std;

//-----------------------------------------------------------------------------------------------
// timediff(double UT1_UTC, double TAI_UTC)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Time differences
 *
 *	@param [in] UT1_UTC
 *	@param [in] TAI_UTC
 *
 *	@return double UT1_TAI
 *	@return double UTC_GPS
 *	@return double UT1_GPS
 *	@return double TT_UTC
 *	@return double GPS_UTC
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<double,double,double,double,double> timediff(double UT1_UTC, double TAI_UTC);

#endif