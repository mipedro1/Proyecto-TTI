// $Source$
//--------------------------------------------------------------------------------
// EccAnom
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file EccAnom.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\EccAnom.hpp"

double EccAnom (double M, double e){
	int maxit = 15;
	int i = 1;
	double eps = DBL_EPSILON;
	// Starting value
	M = fmod(M, 2.0*M_PI);
	
	double E;
	
	if (e<0.8)
		E = M; 
	else
		E = M_PI;
	

	double f = E - e*sin(E) - M;
	E = E - f / ( 1.0 - e*cos(E) );

	// Iteration
	while (fabs(f) > 1e2*eps) {   
		f = E - e*sin(E) - M;
		E = E - f / ( 1.0 - e*cos(E) );
		i = i+1;
		if (i==maxit){
			
			cout<<"convergence problems in EccAnom\n";
			exit(EXIT_FAILURE);
		}  
	}
	return E;

}