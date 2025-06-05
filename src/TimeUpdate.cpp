// $Source$
//--------------------------------------------------------------------------------
// TimeUpdate
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/04
//
/** @file TimeUpdate.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\TimeUpdate.hpp"

Matrix& TimeUpdate(Matrix& P, Matrix& Phi,double Qdt){
	
	Matrix &Phi_t= Phi.transpose();
	P = Phi*P*Phi_t + Qdt;
	return P;

}