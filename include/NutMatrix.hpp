// $Header$
//--------------------------------------------------------------------------------
// NutMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file NutMatrix.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _NutMatrix_
#define _NutMatrix_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"
#include <cmath>

using namespace std;


//-----------------------------------------------------------------------------------------------
// NutMatrix (double Mjd_TT)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Transformation from mean to true equator and equinox
 *
 *	@param [in]  Mjd_TT    Modified Julian Date (Terrestrial Time)
 *
 *	@return Matrix& NutMat    Nutation matrix
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& NutMatrix (double Mjd_TT);

#endif