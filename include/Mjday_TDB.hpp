// $Header$
//--------------------------------------------------------------------------------
// Mjday_TDB
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file Mjday_TDB.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _MjdayTDB_
#define _MjdayTDB_

#include "..\include\matrix.hpp"
#include <cmath>


//-----------------------------------------------------------------------------------------------
// Mjday(int yr, int mon, int day,int hr=0,int min=0,double sec=0)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the Modified Julian Date for barycentric dynamical
 *		   time
 *
 *	@param [in] Mjd_TT      - Modified julian date (TT)
 *
 *	@return double Mjd_TDB     - Modified julian date (TDB)
 *
 */
//-----------------------------------------------------------------------------------------------
double Mjday_TDB(double Mjd_TT);

#endif