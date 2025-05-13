#ifndef _MeasUpdate_
#define _MeasUpdate_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

tuple<Matrix&,Matrix&,Matrix&> MeasUpdate(Matrix& x, double z,double g,double s,Matrix& G,Matrix& P, int n);

#endif