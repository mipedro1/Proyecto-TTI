// $Header$
//--------------------------------------------------------------------------------
// sign_
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file sign_.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Sign_
#define _Sign_


#include <cmath>


//-----------------------------------------------------------------------------------------------
// sign_(double a,double  b)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief returns absolute value of a with sign of b
 *
 *	@param [in] a
 *	@param [in] b
 *
 *	@return double result
 *
 */
//-----------------------------------------------------------------------------------------------
double sign_(double a,double  b);

#endif