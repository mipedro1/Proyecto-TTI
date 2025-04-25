#ifndef _IERS_
#define _IERS_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>

using namespace std;

tuple<double,double,double,double,double,double,double,double,double> IERS(Matrix eop,double Mjd_UTC,char interp);

#endif