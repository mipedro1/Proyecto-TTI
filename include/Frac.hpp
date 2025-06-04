// $Header$
//--------------------------------------------------------------------------------
// Frac
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file Frac.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Frac_
#define _Frac_

#include "..\include\matrix.hpp"
#include <cmath>


//-----------------------------------------------------------------------------------------------
// Frac (double x)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Fractional part of a number (y=x-[x])
 *
 *	@param [in] x
 *
 *	@return double res
 *
 */
//-----------------------------------------------------------------------------------------------
double Frac (double x);

#endif