// $Header$
//--------------------------------------------------------------------------------
// Geodetic
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file Geodetic.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Geodetic_
#define _Geodetic_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Cheb3D.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>


//-----------------------------------------------------------------------------------------------
// Geodetic(Matrix& r)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
 *	       from given position vector (r [m])
 *
 *	@param [in] r         
 *
 *	@return double lon
 *	@return double lat
 *	@return double h 
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<double,double,double> Geodetic(Matrix& r);

#endif