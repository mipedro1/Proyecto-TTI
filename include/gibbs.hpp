// $Header$
//--------------------------------------------------------------------------------
// gibbs
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file gibbs.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Gibbs_
#define _Gibbs_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\angl.hpp"
#include "..\include\unit.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>
#include <string>


//-----------------------------------------------------------------------------------------------
// gibbs( Matrix& r1,Matrix& r2,Matrix& r3)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief this function performs the gibbs method of orbit determination. this
 *	       method determines the velocity at the middle point of the 3 given
 *		   position vectors.
 *
 *	@param [in] r1          - ijk position vector #1         m
 *	@param [in] r2          - ijk position vector #2         m
 *	@param [in] r3          - ijk position vector #3         m
 *
 *	@return Matrix& v2          - ijk velocity vector for r2     m/s
 *	@return double theta       - angl between vectors           rad
 *	@return double theta1
 *	@return double copa
 *	@return String error       - flag indicating success        'ok',...
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<Matrix&,double,double,double,string> gibbs( Matrix& r1,Matrix& r2,Matrix& r3);

#endif