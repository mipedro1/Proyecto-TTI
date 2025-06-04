// $Header$
//--------------------------------------------------------------------------------
// MeasUpdate
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file MeasUpdate.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _MeasUpdate_
#define _MeasUpdate_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;


//-----------------------------------------------------------------------------------------------
// MeasUpdate(Matrix& x, double z,double g,double s,Matrix& G,Matrix& P, int n)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Performs a measurement update step of a Kalman filter.
 *
 *	@param [in,out] x
 *	@param [in] 	z
 *	@param [in] 	g
 *	@param [in] 	s
 *	@param [in] 	G
 *	@param [in,out] P
 *	@param [in] 	n
 *
 *	@return Matrix& K
 *	@return Matrix& x
 *	@return Matrix& P
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix& x, double z,double g,double s,Matrix& G,Matrix& P, int n);

#endif