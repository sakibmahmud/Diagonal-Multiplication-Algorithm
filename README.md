# Diagonal-Multiplication-Algorithm
We present a storage scheme for storing matrices by diagonals and algorithms for performing matrix-matrix and matrix-vector multiplication by diagonals. Matrix elements are accessed with stride-1 and involve no indirect referencing. Access to the transposed matrix requires no additional effort. The proposed storage scheme handles dense matrices and matrices with special structure e.g., banded, triangular, symmetric in a uniform manner. Test results from preliminary numerical experiments with an OpenMP implementation of our method are encouraging.
r8_uniform_01() : this function is used to create inputs. 
diagStore(vector <float> temp) : This function converts and stores a standard matrix (n x n) diagonally to a 1-D vector.
standardMultiply() : Multiplying C= A x B using IJK version. 
standardMultiplyKIJ() : Multiplying C= A x B using KIJ version.
standardMultiplyIKJ() :  Multiplying C= A x B using IKJ version.
diagonalMultiply(int thread, int is_parallel): This function performs the diagonal multiplication C = A x B using diagonal multiplication algorithm. The parallel version of this algorithm uses openmp. You can turn off the parallel code by commenting out pragma block or using is_parallel bool. 
getIndex(): to get start index of a diagonal in a 1-D vector that stores a matrix diagonally. 
