// $Header$
//--------------------------------------------------------------------------------
// hgibbs
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file hgibbs.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Hgibbs_
#define _Hgibbs_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\angl.hpp"
#include "..\include\unit.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>
#include <string>


//-----------------------------------------------------------------------------------------------
// hgibbs( Matrix& r1,Matrix& r2,Matrix& r3,double Mjd1, double Mjd2,double Mjd3)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief this function implements the herrick-gibbs approximation for orbit
 *	       determination, and finds the middle velocity vector for the 3 given
 *		   position vectors.
 *
 *	@param [in] r1          - ijk position vector #1         m
 *	@param [in] r2          - ijk position vector #2         m
 *	@param [in] r3          - ijk position vector #3         m
 *	@param [in] Mjd1        - julian date of 1st sighting    days from 4713 bc
 *	@param [in] Mjd2        - julian date of 2nd sighting    days from 4713 bc
 *	@param [in] Mjd3        - julian date of 3rd sighting    days from 4713 bc
 *
 *	@return Matrix& v2          - ijk velocity vector for r2     m/s
 *	@return double theta       - angl between vectors           rad
 *	@return double theta1
 *	@return double copa
 *	@return String error       - flag indicating success        'ok',...
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<Matrix&,double,double,double,string> hgibbs( Matrix& r1,Matrix& r2,Matrix& r3,double Mjd1, double Mjd2,double Mjd3);

#endif