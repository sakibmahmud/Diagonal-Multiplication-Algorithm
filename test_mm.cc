#include "diagonal_mm.cc"

int main()

{

    
Diagonal_Multiplication m(100);
vector <float> result;
result =m.standardMultiply();
m.standardMultiplyKIJ();
m.diagonalMultiply(2,1);

}
