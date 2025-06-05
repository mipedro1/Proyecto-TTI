// $Header$
//--------------------------------------------------------------------------------
// R_y
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file R_y.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _R_y_
#define _R_y_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// R_y(double angle)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief returns rotation matrix
 *
 *	@param [in] angle       - angle of rotation [rad]
 *
 *	@return Matrix& rotmat      - vector result
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& R_y(double angle);

#endif
