// $Header$
//--------------------------------------------------------------------------------
// global
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file global.hpp
 *  @brief This header file contains auxiliar data used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#ifndef _GLOBAL_
#define _GLOBAL_

#include "../include/matrix.hpp"
#include "../include/Mjday.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <string.h>


typedef struct{
	double Mjd_UTC,Mjd_TT;
	int n,m,sun,moon,planets;
} Param;

extern Param AuxParam;
extern Matrix eopdata;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix PC;
extern Matrix obs;

void eop19620101(int c);
void GGM03S(int n);
void DE430Coeff(int f, int c);
void GEOS3(int f);

#endif