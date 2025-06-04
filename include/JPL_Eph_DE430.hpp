// $Header$
//--------------------------------------------------------------------------------
// JPL_Eph_DE430
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file JPL_Eph_DE430.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _JplEp_
#define _JplEp_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Cheb3D.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;


//-----------------------------------------------------------------------------------------------
// JPL_Eph_DE430(double Mjd_TDB)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the sun, moon, and nine major planets' equatorial
 *	       position using JPL Ephemerides
 *
 *	@param [in] Mjd_TDB         Modified julian date of TDB
 *
 *	@return Matrix& r_Earth(solar system barycenter (SSB))
 *	@return Matrix& r_Mars
 *	@return Matrix& r_Mercury
 *	@return Matrix& r_Venus
 *	@return Matrix& r_Jupiter
 *	@return Matrix& r_Saturn	
 *	@return Matrix& r_Uranus	
 *	@return Matrix& r_Neptune
 *	@return Matrix& r_Pluto	
 *	@return Matrix& r_Moon
 *	@return Matrix& r_Sun(geocentric equatorial position ([m]) referred to the
 *                  International Celestial Reference Frame (ICRF))
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&> JPL_Eph_DE430(double Mjd_TDB);

#endif