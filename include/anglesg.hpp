// $Header$
//--------------------------------------------------------------------------------
// anglesg
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file anglesg.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Anglesg_
#define _Anglesg_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\angl.hpp"
#include "..\include\unit.hpp"
#include "..\include\Geodetic.hpp"
#include "..\include\LTC.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\gibbs.hpp"
#include "..\include\hgibbs.hpp"
#include "..\include\elements.hpp"
#include "..\include\rpoly.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>
#include <string>

using namespace std;

//-----------------------------------------------------------------------------------------------
// anglesg ( double az1,double az2,double az3,double el1,double el2,double el3,double Mjd1,double Mjd2,double Mjd3,Matrix& Rs1,Matrix& Rs2,Matrix& Rs3 )
//-----------------------------------------------------------------------------------------------
/**
 *	@brief this function solves the problem of orbit determination using three
 *	       optical sightings.
 *
 *	@param [in] az1      - azimuth at t1               rad
 *	@param [in] az2      - azimuth at t2               rad
 *	@param [in] az3      - azimuth at t3               rad
 *	@param [in] el1      - elevation at t1             rad
 *	@param [in] el2      - elevation at t2             rad
 *	@param [in] el3      - elevation at t3             rad
 *	@param [in] Mjd1     - Modified julian date of t1
 *	@param [in] Mjd2     - Modified julian date of t2
 *	@param [in] Mjd3     - Modified julian date of t3
 *	@param [in] Rs1      - ijk site1 position vector   m
 *	@param [in] Rs2      - ijk site2 position vector   m
 *	@param [in] Rs3      - ijk site3 position vector   m
 *
 *	@return Matrix& r        - ijk position vector at t2   m
 *	@return Matrix& v        - ijk velocity vector at t2   m/s
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<Matrix&,Matrix&> anglesg ( double az1,double az2,double az3,double el1,double el2,double el3,double Mjd1,double Mjd2,double Mjd3,Matrix& Rs1,Matrix& Rs2,Matrix& Rs3 );

#endif