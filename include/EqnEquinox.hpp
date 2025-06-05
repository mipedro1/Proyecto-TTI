// $Header$
//--------------------------------------------------------------------------------
// EqnEquinox
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file EqnEquinox.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _EqnEquinox_
#define _EqnEquinox_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include <cmath>

using namespace std;

//-----------------------------------------------------------------------------------------------
// EqnEquinox (double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computation of the equation of the equinoxes
 *
 *	@param [in] Mjd_TT    Modified Julian Date (Terrestrial Time)   
 *
 *	@return double EqE      Equation of the equinoxes
 *
 */
//-----------------------------------------------------------------------------------------------
double EqnEquinox (double Mjd_TT);

#endif