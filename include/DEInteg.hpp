// $Header$
//--------------------------------------------------------------------------------
// DEInteg
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file DEInteg.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _DEInteg_
#define _DEInteg_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\sign_.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

//-----------------------------------------------------------------------------------------------
// DEInteg(Matrix& f (double t,Matrix& y),double t,double tout,double relerr,double abserr,int n_eqn,Matrix& y)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Chebyshev approximation of 3-dimensional vectors
 *
 *	@param [in] f       
 *	@param [in] t      
 *	@param [in] tout      
 *	@param [in] relerr      
 *	@param [in] abserr      
 *	@param [in] n_eqn      
 *	@param [in,out] y      
 *
 *	@return Matrix& y   
 *
 */
//-----------------------------------------------------------------------------------------------
 Matrix& DEInteg(Matrix& f (double t,Matrix& y),double t,double tout,double relerr,double abserr,int n_eqn,Matrix& y);

#endif