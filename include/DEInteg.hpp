#ifndef _DEInteg_
#define _DEInteg_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\IERS.hpp"
#include "..\include\timediff.hpp"
#include "..\include\PrecMatrix.hpp"
#include "..\include\NutMatrix.hpp"
#include "..\include\PoleMatrix.hpp"
#include "..\include\GHAMatrix.hpp"
#include "..\include\JPL_Eph_DE430.hpp"
#include "..\include\AccelHarmonic.hpp"
#include "..\include\G_AccelHarmonic.hpp"
#include "..\include\AccelPointMass.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "..\include\sign_.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

 Matrix& DEInteg(Matrix& f (double t,Matrix& y),double t,double tout,double relerr,double abserr,int n_eqn,Matrix& y);

#endif