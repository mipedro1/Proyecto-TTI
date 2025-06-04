// $Header$
//--------------------------------------------------------------------------------
// Accel
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file Accel.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Accel_
#define _Accel_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

//-----------------------------------------------------------------------------------------------
// Accel(double x,Matrix& Y)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the acceleration of an Earth orbiting satellite due to 
 *         - the Earth's harmonic gravity field, 
 *         - the gravitational perturbations of the Sun and Moon,
 *         - the solar radiation pressure, and
 *         - the atmospheric drag.
 *
 *	@param [in] x    Terrestrial Time (Modified Julian Date).
 *	@param [in] Y Satellite state vector in the ICRF/EME2000 system.
 *
 *	@return Matrix& dY  Acceleration (a = d^2r/dt^2) in the ICRF/EME2000 system.
 *
 */
//-----------------------------------------------------------------------------------------------
 Matrix& Accel(double x,Matrix& Y);

#endif