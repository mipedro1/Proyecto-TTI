// $Header$
//--------------------------------------------------------------------------------
// G_AccelHarmonic
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file G_AccelHarmonic.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _GAccel_
#define _GAccel_


#include "..\include\matrix.hpp"
#include "..\include\AccelHarmonic.hpp"
#include <cmath>

using namespace std;


//-----------------------------------------------------------------------------------------------
// G_AccelHarmonic( Matrix& r,Matrix& U, int n_max,int  m_max )
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the gradient of the Earth's harmonic gravity field 
 *
 *	@param [in] r           Satellite position vector in the true-of-date system 
 *	@param [in] U           Transformation matrix to body-fixed system
 *	@param [in] n           Gravity model degree
 *	@param [in] m 			Gravity model order
 *
 *	@return Matrix& G    		Gradient (G=da/dr) in the true-of-date system
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& G_AccelHarmonic( Matrix& r,Matrix& U, int n_max,int  m_max );

#endif