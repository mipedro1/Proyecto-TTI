#ifndef _JplEp_
#define _JplEp_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "..\include\Cheb3D.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

tuple<Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&,Matrix&> JPL_Eph_DE430(double Mjd_TDB);

#endif