
#include "diagonal_mv.h"

Matrix_Vector::Matrix_Vector(int n)
{
	dim=n;
	int seed=3145;
	cout << dim <<endl;
     A.resize(dim*dim);
	// taking input to Matrix A
        for (int i = 0; i < dim; i++) {

                for (int j = 0; j < dim; j++) {

                        A[i * n + j] =r8_uniform_01( & seed);
                        
                }

        }

    //taking input into vector X
    
    int seed1=101112131;
    X.resize(dim);

     for(int i=0; i <dim; i++)
     {
                      
            X[i]= r8_uniform_01(&seed1);
                           
     }   

     //Storing Matrix A in diagonal order

      diagMatrixA = diagStore(A);

      //Storing the index order of diagonal
      //For example, index of main diagonal is 0, index of first subdiagonal -1
       storeIndex();
 
	
}

void Matrix_Vector::storeIndex()
{
	//atmost 2n-1 diagonals in a n x n matrix
	diagonalIndex.resize(2*dim-1);
	 for(int i=0; i<dim; i++){

           diagonalIndex[i]=i;

        }
   int k=1;
   for(int i=dim; i< 2*dim-1; i++){

          diagonalIndex[i]=k*(-1);
          k++;
      } 



}

//Performing y=Ax in terms of diagonals

vector<float> Matrix_Vector::matrixvectorMultiply()
{ 
    vector<float> Y;
    Y.resize(dim);
    int n=dim;

 //Computing with Main and Superdiagonal
   for(int d=0; d<2*dim-1; d++)
   { //loop 1 begins
          
        int k=diagonalIndex[d];
        int temp=0;
      
          if(k>=0)
          	{ 
   
              int i=k;
              int j=0;
     
               int start_index=k*n-k*(k-1)/2;

               int stop_index=start_index+(n-k)-1;
               

             for(; start_index<=stop_index;start_index++)
             { //loop 2 begins

                    Y[j]+=A[start_index]*X[i];
                        
                   
                       j++;

                       i++;

              } //loop 2 ends

           } 
           

//Computing Sub diagonals
         else

          {

                int abs_k = abs(k);

               int i=0;
               int j=abs_k;

               int start_index=n*(n+1)/2 + (abs_k-1)*n-abs_k*(abs_k-1)/2;
               int stop_index=start_index+n-abs_k-1;
               for(;start_index<=stop_index; start_index++)
               	{ //loop 3 begins

                         Y[j]+=A[start_index]*X[i];

                          j++;

                          i++;

                 } //loop 3 ends
              

           }
         
      } //loop 1 ends

  return Y;


}

//Performa y= A^Tx in terms of diagonals


vector<float> Matrix_Vector::matrixvectorTransposeMultiply()
{

  vector<float> Y;
    Y.resize(dim);
    int n=dim;

 //Computing with Main and Superdiagonal
   for(int d=0; d<2*dim-1; d++)
   { //loop 1 begins
          
        int k=diagonalIndex[d];
        int temp=0;
      
          if(k>=0)
            { 
             //superdiagonal become subdiagonals and vice verca 
              int i=0;
              int j=k;
     
               int start_index=k*n-k*(k-1)/2;

               int stop_index=start_index+(n-k)-1;
               

             for(; start_index<=stop_index;start_index++)
             { //loop 2 begins

                    Y[j]+=A[start_index]*X[i];
                        
                   
                       j++;

                       i++;

              } //loop 2 ends

           } 
           

//Computing Sub diagonals
         else

          {

                int abs_k = abs(k);

               int i=abs_k;
               int j=0;

               int start_index=n*(n+1)/2 + (abs_k-1)*n-abs_k*(abs_k-1)/2;
               int stop_index=start_index+n-abs_k-1;
               for(;start_index<=stop_index; start_index++)
                { //loop 3 begins

                         Y[j]+=A[start_index]*X[i];

                          j++;

                          i++;

                 } //loop 3 ends
              

           }
         
      } //loop 1 ends

  return Y;






}




