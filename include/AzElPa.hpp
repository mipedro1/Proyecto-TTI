// $Header$
//--------------------------------------------------------------------------------
// AzElPa
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file AzElPa.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _AzElPa_
#define _AzElPa_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>

using namespace std;

//-----------------------------------------------------------------------------------------------
// AzElPa(Matrix& s)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes azimuth, elevation and partials from local tangent coordinates
 *
 *	@param [in] s      Topocentric local tangent coordinates (East-North-Zenith frame)
 *
 *	@return double A      Azimuth [rad]
 *	@return double E      Elevation [rad]
 *	@return Matrix& dAds   Partials of azimuth w.r.t. s
 *	@return Matrix& dEds   Partials of elevation w.r.t. s
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix& s);

#endif