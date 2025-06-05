// $Header$
//--------------------------------------------------------------------------------
// PrecMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file PrecMatrix.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _PrecMatrix_
#define _PrecMatrix_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_z.hpp"
#include <cmath>

using namespace std;


//-----------------------------------------------------------------------------------------------
// PrecMatrix (double Mjd_1,double  Mjd_2)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Precession transformation of equatorial coordinates
 *
 *	@param [in] Mjd_1     Epoch given (Modified Julian Date TT)
 *	@param [in] MjD_2     Epoch to precess to (Modified Julian Date TT)
 *
 *	@return Matrix& PrecMat   Precession transformation matrix
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& PrecMatrix (double Mjd_1,double  Mjd_2);

#endif