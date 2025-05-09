#include "..\include\TimeUpdate.hpp"

Matrix& TimeUpdate(Matrix& P, Matrix& Phi,double Qdt){
	
	Matrix &Phi_t= Phi.transpose();
	P = Phi*P*Phi_t + Qdt;
	return P;

}