// $Source$
//--------------------------------------------------------------------------------
// angl
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file angl.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//-------------------------------------------------------------------------------- 
#include "..\include\angl.hpp"

double angl ( Matrix& vec1, Matrix& vec2 ){
	double small     = 0.00000001;
	double undefined = 999999.1;

	double magv1 = norm(vec1);
	double magv2 = norm(vec2);
	
	double eps=2.22044604925031e-16;
	
	double theta;
	if (magv1*magv2 > pow(small,2)){
		double temp= vec1.dot(vec2) / (magv1*magv2);
		if (fabs( temp ) > 1.0)
			temp= sign_(temp, temp) * 1.0;
		
		theta= acos( temp );
	}else
		theta= eps;
	
	
	return theta;
}