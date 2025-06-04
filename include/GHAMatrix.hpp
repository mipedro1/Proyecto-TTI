// $Header$
//--------------------------------------------------------------------------------
// GHAMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file GHAMatrix.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _GHAMatrix_
#define _GHAMatrix_


#include "..\include\matrix.hpp"
#include "..\include\gast.hpp"
#include "..\include\R_z.hpp"
#include <cmath>

using namespace std;


//-----------------------------------------------------------------------------------------------
// GHAMatrix (double Mjd_UT1)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Transformation from true equator and equinox to Earth equator and 
 *		   Greenwich meridian system 
 *
 *	@param [in] Mjd_UT1   Modified Julian Date UT1
 *
 *	@return Matrix& GHAmat    Greenwich Hour Angle matrix
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& GHAMatrix (double Mjd_UT1);

#endif