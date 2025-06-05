// $Source$
//--------------------------------------------------------------------------------
// gibbs
//--------------------------------------------------------------------------------
// Proyecto-TTI
//
// Copyright (c) 2020, Meysam Mahooti
//
// Created: 2025/06/05
//
/** @file gibbs.cpp
 *  @brief Auxiliar function used by EKF_GEOS3
 *
 *	@author Miguel de Pedro Olagaray
 *	@bug No knows bugs.
 */
//--------------------------------------------------------------------------------
#include "..\include\gibbs.hpp"

tuple<Matrix&,double,double,double,string> gibbs( Matrix& r1,Matrix& r2,Matrix& r3){
	double small= 0.00000001;
	double theta= 0.0;
	string error = "          ok";
	double theta1= 0.0;

	double magr1 = norm( r1 );
	double magr2 = norm( r2 );
	double magr3 = norm( r3 );
	Matrix& v2=zeros(3);
	for (int i=1;i<=3;i++){
		v2(i)= 0.0;
	}

	Matrix& p = r2.cross( r3 );
	Matrix& q = r3.cross( r1 );
	Matrix& w = r1.cross( r2 );
	Matrix& pn = unit( p );
	Matrix& r1n = unit( r1 );
	double copa=  asin( pn.dot( r1n ) );

	if ( fabs( r1n.dot(pn) ) > 0.017452406 )  
		error= "not coplanar";
	

	Matrix& d = p + q + w;
	double magd = norm(d);
	Matrix& n = p*magr1 + q*magr2 + w*magr3;
	double magn = norm(n);
	Matrix& nn = unit( n );
	Matrix& dn = unit( d );

	// -------------------------------------------------------------
	// determine if  the orbit is possible. both d and n must be in
	// the same direction, and non-zero.
	// -------------------------------------------------------------
	if ( ( fabs(magd)<small ) || ( fabs(magn)<small ) || 
	   ( nn.dot(dn) < small ) ){
		error= "  impossible";
	   }else{
		  theta  = angl( r1,r2 );
		  theta1 = angl( r2,r3 );

		  // ----------- perform gibbs method to find v2 -----------
		  double r1mr2= magr1-magr2;
		  double r3mr1= magr3-magr1;
		  double r2mr3= magr2-magr3;
		  Matrix& s  = r3*r1mr2 + r2*r3mr1 + r1*r2mr3;
		  Matrix& b  = d.cross( r2 );
		  double l  = sqrt(GM_Earth / (magd*magn) );
		  double tover2= l / magr2;
		  v2 = b*tover2  + s*l ;
    }
	return tie(v2, theta,theta1,copa, error);
}