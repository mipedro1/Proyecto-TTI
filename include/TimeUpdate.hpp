#ifndef _TimeUpdate_
#define _TimeUpdate_


#include "..\include\matrix.hpp"
#include <cmath>

using namespace std;

Matrix& TimeUpdate(Matrix& P, Matrix& Phi,double Qdt=0.0);

#endif