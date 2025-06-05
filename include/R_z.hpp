// $Header$
//--------------------------------------------------------------------------------
// R_z
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file R_z.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _R_z_
#define _R_z_

#include "..\include\matrix.hpp"
#include <cmath>

//-----------------------------------------------------------------------------------------------
// R_z(double angle)
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
Matrix& R_z(double angle);

#endif
