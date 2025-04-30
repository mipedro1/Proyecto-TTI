#ifndef _AzElPa_
#define _AzElPa_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include <cmath>
#include <tuple>

using namespace std;

tuple<double,double,Matrix&,Matrix&> AzElPa(Matrix& s);

#endif