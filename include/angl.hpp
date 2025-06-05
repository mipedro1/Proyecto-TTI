// $Header$
//--------------------------------------------------------------------------------
// angl
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file angl.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Angl_
#define _Angl_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\sign_.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>


//-----------------------------------------------------------------------------------------------
// angl ( Matrix& vec1, Matrix& vec2 )
//-----------------------------------------------------------------------------------------------
/**
 *	@brief this function returns theta, an angle between the two vectors
 *
 *	@param [in] vec1         - vector 1   
 *	@param [in] vec2         - vector 2
 *
 *	@return double theta        - angle between the two vectors  -pi to pi
 *
 */
//-----------------------------------------------------------------------------------------------
double angl ( Matrix& vec1, Matrix& vec2 );

#endif