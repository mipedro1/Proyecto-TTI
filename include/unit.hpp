// $Header$
//--------------------------------------------------------------------------------
// unit
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file unit.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Unit_
#define _Unit_


#include "..\include\matrix.hpp"
#include <cmath>
#include <tuple>


//-----------------------------------------------------------------------------------------------
// unit(Matrix& vec)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief this function calculates a unit vector given the original vector. if a
 *	       zero vector is input, the vector is set to zero.
 *
 *	@param [in] vec         - vector
 *
 *	@return Matrix& outvec      - unit vector
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& unit(Matrix& vec);

#endif