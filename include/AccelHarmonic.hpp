// $Header$
//--------------------------------------------------------------------------------
// AccelHarmonic
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file AccelHarmonic.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _AccelHarmonic_
#define _AccelHarmonic_


#include "..\include\matrix.hpp"
#include "../include/global.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Legendre.hpp"
#include <cmath>

using namespace std;
//-----------------------------------------------------------------------------------------------
// AccelHarmonic(Matrix& r, Matrix& E, double n_max,double m_max)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the acceleration due to the harmonic gravity field of the central body.
 *
 *	@param [in] r           Satellite position vector in the inertial system
 *	@param [in] E           Transformation matrix to body-fixed system
 *	@param [in] n_max       Maximum degree
 *	@param [in] m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
 *
 *	@return Matrix& a           Acceleration (a=d^2r/dt^2)
 *
 */
//-----------------------------------------------------------------------------------------------
 Matrix& AccelHarmonic(Matrix& r, Matrix& E, double n_max,double m_max);

#endif