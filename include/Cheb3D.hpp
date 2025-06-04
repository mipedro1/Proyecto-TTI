// $Header$
//--------------------------------------------------------------------------------
// Cheb3D
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file Cheb3D.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Cheb3D_
#define _Cheb3D_

#include "..\include\matrix.hpp"
#include <cmath>


//-----------------------------------------------------------------------------------------------
// Cheb3D(double t, int N, double Ta,double Tb, Matrix& Cx, Matrix& Cy,Matrix& Cz)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Chebyshev approximation of 3-dimensional vectors
 *
 *	@param [in] N       Number of coefficients
 *	@param [in] Ta      Begin interval
 *	@param [in] Tb      End interval
 *	@param [in] Cx      Coefficients of Chebyshev polyomial (x-coordinate)
 *	@param [in] Cy      Coefficients of Chebyshev polyomial (y-coordinate)
 *	@param [in] Cz      Coefficients of Chebyshev polyomial (z-coordinate)
 *
 *	@return Matrix& ChebApp   Chebyshev approximation of 3-dimensional vectors
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& Cheb3D(double t, int N, double Ta,double Tb, Matrix& Cx, Matrix& Cy,Matrix& Cz);

#endif