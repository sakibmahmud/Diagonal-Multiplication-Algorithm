
#ifndef DIAGONAL_H
#define DIAGONAL_H
#include <iostream>
#include <vector>
#include <omp.h>
#include <cmath> 
#include <ctime> 
#include <cstdlib> 
#include <iterator> 
#include <algorithm> 
#include <csignal> 
using namespace std;

class Diagonal_Multiplication
{
 public:
   Diagonal_Multiplication(int);

   vector<float> standardMultiply();
   vector<float> standardMultiplyKIJ();
   vector<float> standardMultiplyIKJ();
   vector<float> diagonalMultiply(int, int);
   float r8_uniform_01(int*);
   vector<float> diagStore(vector<float>);
   int getIndex(int);

 private:
   int dim;
   vector <float> A;
   vector <float> B;
   vector<float> diagMatrixA;
   vector<float> diagMatrixB;
   


};

#endif