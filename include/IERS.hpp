// $Header$
//--------------------------------------------------------------------------------
// IERS
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file IERS.h
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _IERS_
#define _IERS_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

//-----------------------------------------------------------------------------------------------
// IERS(double Mjd_UTC,char interp='n')
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Management of IERS time and polar motion data
 *
 *	@param [in] Mjd_UTC
 *	@param [in] interp
 *
 *	@return double x_pole
 *	@return double y_pole    		
 *	@return double UT1_UTC    		
 *	@return double LOD    		
 *	@return double dpsi    		
 *	@return double deps    		
 *	@return double dx_pole    		
 *	@return double dy_pole    		
 *	@return double TAI_UTC    		
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<double,double,double,double,double,double,double,double,double> IERS(double Mjd_UTC,char interp='n');

#endif