// $Header$
//--------------------------------------------------------------------------------
// Mjday
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file Mjday.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Mjday_
#define _Mjday_

#include "..\include\matrix.hpp"
#include <cmath>


//-----------------------------------------------------------------------------------------------
// Mjday(int yr, int mon, int day,int hr=0,int min=0,double sec=0)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Converts a given date and time to Modified Julian Date (MJD)
 *
 *	@param [in] year        - year
 *	@param [in] mon         - month
 *	@param [in] day         - day 
 *	@param [in] hr          - universal time hour 
 *	@param [in] min         - universal time min 
 *	@param [in] sec         - universal time sec 
 *
 *	@return double Mjd         - Modified julian date      
 *
 */
//-----------------------------------------------------------------------------------------------
double Mjday(int yr, int mon, int day,int hr=0,int min=0,double sec=0);

#endif