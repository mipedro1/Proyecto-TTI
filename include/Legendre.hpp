// $Header$
//--------------------------------------------------------------------------------
// Legendre
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file Legendre.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Legendre_
#define _Legendre_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;



//-----------------------------------------------------------------------------------------------
// Legendre(int n,int m,double fi)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the Legendre polynomials and their derivatives.
 *
 *	@param [in] n
 *	@param [in] m
 *	@param [in] fi
 *
 *	@return Matrix& pnm
 *	@return Matrix& dpnm
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<Matrix&,Matrix&> Legendre(int n,int m,double fi);

#endif