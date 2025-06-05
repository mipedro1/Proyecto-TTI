// $Header$
//--------------------------------------------------------------------------------
// VarEqn
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file VarEqn.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _VarEqn_
#define _VarEqn_


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
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

//-----------------------------------------------------------------------------------------------
// VarEqn(double x,Matrix& yPhi)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the variational equations, i.e. the derivative of the state vector
 *	       and the state transition matrix
 *
 *	@param [in] x           Time since epoch in [s]
 *	@param [in] yPhi        (6+36)-dim vector comprising the state vector (y) and the
 *							state transition matrix (Phi) in column wise storage order
 *
 *	@return Matrix& yPhip       Derivative of yPhi
 *
 */
//-----------------------------------------------------------------------------------------------
 Matrix& VarEqn(double x,Matrix& yPhi);

#endif