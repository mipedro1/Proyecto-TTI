// $Source$
//--------------------------------------------------------------------------------
// unit
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file unit.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//-------------------------------------------------------------------------------- 
#include "..\include\unit.hpp"

Matrix& unit(Matrix& vec){
	
	double small = 0.000001;
	double magv = norm(vec);
	
	Matrix& outvec=zeros(3);
	
	if ( magv > small )
		for (int i=1;i<=3;i++){
			outvec(i)= vec(i)/magv;
		}
	else
		for (int i=1;i<=3;i++){
			outvec(i)= 0.0;
		}
	
	return outvec;
}
