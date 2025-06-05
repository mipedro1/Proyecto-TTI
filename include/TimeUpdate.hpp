// $Header$
//--------------------------------------------------------------------------------
// TimeUpdate
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file TimeUpdate.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _TimeUpdate_
#define _TimeUpdate_


#include "..\include\matrix.hpp"
#include <cmath>

using namespace std;

//-----------------------------------------------------------------------------------------------
// TimeUpdate(Matrix& P, Matrix& Phi,double Qdt=0.0)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief time update
 *
 *	@param [in] P
 *	@param [in] Phi
 *	@param [in] Qdt
 *
 *	@return Matrix& P
 *
 */
//-----------------------------------------------------------------------------------------------
Matrix& TimeUpdate(Matrix& P, Matrix& Phi,double Qdt=0.0);

#endif