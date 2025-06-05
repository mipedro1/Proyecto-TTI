// $Header$
//--------------------------------------------------------------------------------
// MeanObliquity
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file MeanObliquity.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _MeanObliquity_
#define _MeanObliquity_

#include "..\include\matrix.hpp"
#include <cmath>
#include "..\include\SAT_Const.hpp"


//-----------------------------------------------------------------------------------------------
// MeanObliquity (double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the mean obliquity of the ecliptic
 *
 *	@param [in] Mjd_TT    Modified Julian Date (Terrestrial Time)
 *
 *	@return double MOblq     Mean obliquity of the ecliptic [rad]
 *
 */
//-----------------------------------------------------------------------------------------------
double MeanObliquity (double Mjd_TT);

#endif