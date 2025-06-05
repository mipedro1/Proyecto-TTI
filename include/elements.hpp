// $Header$
//--------------------------------------------------------------------------------
// elements
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file elements.hpp
 *  @brief This header file contains an auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _Elements_
#define _Elements_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Cheb3D.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>


//-----------------------------------------------------------------------------------------------
// elements (Matrix& y)
//-----------------------------------------------------------------------------------------------
/**
 *	@brief Computes the osculating Keplerian elements from the satellite state
 *	       vector for elliptic orbits
 *
 *	@param [in] y        State vector (x,y,z,vx,vy,vz)   
 *
 *	@return double p        semilatus rectum [m]
 *	@return double a        Semimajor axis 
 *	@return double e        Eccentricity 
 *	@return double i        Inclination [rad]
 *	@return double Omega    Longitude of the ascending node [rad]
 *	@return double omega    Argument of pericenter [rad]
 *	@return double M        Mean anomaly [rad]
 *
 */
//-----------------------------------------------------------------------------------------------
tuple<double,double,double,double,double,double,double> elements (Matrix& y);

#endif