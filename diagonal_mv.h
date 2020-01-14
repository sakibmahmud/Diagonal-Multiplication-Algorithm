#ifndef MVECTOR_H
#define MVECTOR_H

#include <iostream>
#include <cmath>
#include "diagonal_mm.cc"
using namespace std;

//Derive Class: Base Class
class Matrix_Vector:public Diagonal_Multiplication
{
public:
	Matrix_Vector(int);
	vector<float> matrixvectorMultiply();
	vector<float> matrixvectorTransposeMultiply();
	void storeIndex();

private:
	vector<float> X;
	vector<int> diagonalIndex;
  
	

};

#endif