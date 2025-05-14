#ifndef _Accel_
#define _Accel_


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
#include "..\include\AccelPointMass.hpp"
#include "..\include\Mjday_TDB.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

 Matrix& Accel(double x,Matrix& Y);

#endif