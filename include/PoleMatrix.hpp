// $Header$
//--------------------------------------------------------------------------------
// PoleMatrix
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file PoleMatrix.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _PoleMatrix_
#define _PoleMatrix_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include "..\include\R_y.hpp"
#include "..\include\R_x.hpp"
#include <cmath>

using namespace std;


//-----------------------------------------------------------------------------------------------
// PoleMatrix (double xp,double yp)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates
 *		   for a given date
 *
 *	@param [in]  xp
 *	@param [in]  yp   
 *  Pole coordinte(xp,yp)			 
 *
 *	@return Matrix& PoleMat   Pole matrix
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& PoleMatrix (double xp,double yp);

#endif