#ifndef _NutMatrix_
#define _NutMatrix_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include "..\include\MeanObliquity.hpp"
#include "..\include\NutAngles.hpp"
#include "..\include\R_x.hpp"
#include "..\include\R_z.hpp"
#include <cmath>

using namespace std;

Matrix& NutMatrix (double Mjd_TT);

#endif