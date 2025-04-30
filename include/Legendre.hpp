#ifndef _Legendre_
#define _Legendre_


#include "..\include\matrix.hpp"
#include "..\include\SAT_Const.hpp"
#include "../include/global.hpp"
#include <cmath>
#include <tuple>

using namespace std;

tuple<Matrix&,Matrix&> Legendre(int n,int m,double fi);

#endif