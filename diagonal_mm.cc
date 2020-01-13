#include "diagonal_mm.h"

//generating random values for matrix input
float Diagonal_Multiplication::r8_uniform_01(int *seed)
{

 int k;
		float r;

		k = * seed / 127773;

		* seed = 16807 * ( * seed - k * 127773) - k * 2836;

		if ( * seed < 0) {
		        * seed = * seed + 2147483647;
		}

		r = (float)( * seed) * 4.656612875E-10;

		return r;


}

//Constructor
Diagonal_Multiplication::Diagonal_Multiplication(int n)
{

	cout << n <<endl;
	n=dim;
    int seed = 32456;
    A.resize(n*n);
    B.resize(n*n);
	// taking input to A
        for (int i = 0; i < dim; i++) {

                for (int j = 0; j < dim; j++) {

                        A[i * n + j] =r8_uniform_01( & seed);
                        
                }

        }

        //storing matrix A diagonally in diagMatrixA;
        diagMatrixA = diagStore(A);

      //taking input into B
      
      for (int i = 0; i < n; i++) {

                for (int j = 0; j < n; j++) {

                       B[i * n + j] = r8_uniform_01( & seed);
                       

                }

        }  

        //storing matrix B diagonally in diagMatrixB
        diagMatrixB=diagStore(B);
}

//storing a n x n matrix in diagonal order: Main-Super-Sub diagonal
vector<float> Diagonal_Multiplication::diagStore(vector <float> temp)
{
   vector<float> diag_st;
   diag_st.resize(dim*dim);

   int index1 = 0;

        for (int i = 0; i < dim; i++) { //improve here length1=n-1 since <length1, length1 is equal to total number of main and super diagonals here which is k
                for (int j = 0; j < dim - i; j++) { //improve here //improve here this is mth component for ak(m)=a_m,k+m = a[m][k+m]

                        diag_st[index1] = temp[j * dim + (j + i)]; //j*n+(j+i)

                        index1++;

                }

        }

        int index2 = dim * (dim + 1) / 2;

        // cout << index2 <<endl;

        for (int i = 0; i < dim - 1; i++) {
                for (int j = 1; j < dim - i; j++) {

                        diag_st[index2] = temp[(j + i) * dim + (j - 1)]; //(j+i)*n+(j-1)

                        index2++;

                }

        }

        return diag_st;

}

//Multiplying C= A X B and returning results 
//IJK version
vector<float> Diagonal_Multiplication::standardMultiply()
{
  vector <float> C;
  C.resize(dim*dim);

    for (int i = 0; i < dim; i++)
               {
             for (int j = 0; j < dim; j++ )
    		  {
     
                    //time = clock();
                for (int k = 0; k < dim; k++ )
                      {
                       

                          C[i*dim+j]=C[i*dim+j]+A[i*dim+k]*B[k*dim+j];
                      }

                   
                 }
              }

  return C;
}

//KIJ version
vector <float> Diagonal_Multiplication::standardMultiplyKIJ()
{
	vector <float> C;
	C.resize(dim*dim);

	 for (int k = 0; k < dim; k++)
      {
             for (int i = 0; i < dim; i++ )
    		  {
     
                  
                for (int j = 0; j < dim; j++ )
                      {
                       

                          C[i*dim+j]=C[i*dim+j]+A[i*dim+k]*B[k*dim+j];
                      }

                   

                 }
              }

             return C; 
}

//IKJ version

vector<float> Diagonal_Multiplication::standardMultiplyIKJ()
{

   vector<float> C;
   C.resize(dim*dim);

    for (int i = 0; i < dim; i++)
               
        {
             for (int k = 0; k < dim; k++ )
    		  {
    
                for (int j = 0; j < dim; j++ )
                      {
                       

                          C[i*dim+j]=C[i*dim+j]+A[i*dim+k]*B[k*dim+j];
                      }

                  

                 }
              }

   return C;
}

//C=AB using diagonal multiplication algorithm
vector<float> Diagonal_Multiplication::diagonalMultiply(int thread, int is_parallel)
{

  vector<float> C;
  int i,j,k;
  int n=dim;
  C.resize(dim*dim);
  float test=0.0; 
 // if(is_parallel==1)
  omp_set_num_threads(thread);

  //computing super-diagonal and main-diagonal

  #pragma omp parallel shared(C, diagMatrixA, diagMatrixB) private(i,j,k)
        {
          
           
           int start_index1, start_index2, start_index3, start_index4, start_index5, start_index6, start_index;
          
          #pragma omp for schedule(dynamic) reduction(+:test)
        for (k = 0; k <= n - 1; k++) { //1 

                start_index = k * n - k * (k - 1) / 2;
                
              
                for ( i = k + 1; i <= n - 1; i++) { //2

                       
                         start_index2 = getIndex(k-i);

                         start_index1 = getIndex(i);
                        
                        int length1 = n - i;

                       // float start1 = clock();
                        for (j = 0; j <= length1 - 1; j++) { //3 
                                //cout << "OK";
                                 test =  C[start_index + i - k + j] + diagMatrixA[start_index2] * diagMatrixB[start_index1];
                                  C[start_index + i - k + j] =test;
                                start_index1++;
                                start_index2++;

                        } //3
                         // float end1 = clock();
                       // total_time = total_time + (end1 - start1);
                   
                         start_index3 = getIndex(i);
                       
                        

                         start_index4 = getIndex(k-i);
                       
                

                       // float start2 = clock();
                        for (j = 0; j <= length1 - 1; j++) { //4
                              
                                 test =C[start_index + j] +  diagMatrixA[start_index3] * diagMatrixB[start_index4 + k];

                                C[start_index + j] = test;
                 
                                start_index3++;
                                start_index4++;

                        } //4

                      //  float end2 = clock();

                      //  total_time = total_time + (end2 - start2);

                } //2

                for ( i = 0; i <= k; i++) { //5
 
                      
                         start_index5 = getIndex(i);
                        
                         start_index6 =getIndex(k-i);

                        int length3 = n - i;
                        length3 = length3 - (k - i);

                      //  float start3 = clock();

                        for (j = 0; j <= length3 - 1; j++) { //6

                                test = C[start_index + j] + diagMatrixA[start_index5] * diagMatrixB[start_index6 + i];

                                C[start_index + j] = test;
                                start_index5++;
                                start_index6++;

                                
                        } //6

                        //float end3 = clock();
                       // total_time = total_time + (end3 - start3);

                } //5

        } //1

 } //parallel block

        // computing sub diagonal
      #pragma omp parallel shared(C, diagMatrixA, diagMatrixB) private(i,j,k)
      {
           
           int start_index7, start_index8, start_index9, start_index10, start_index11, start_index12, start_index;
          
         #pragma omp for schedule(dynamic) reduction(+:test)
        for ( k = 1; k <= n - 1; k++) { //7

                start_index = n * (n + 1) / 2 + (k - 1) * n - k * (k - 1) / 2;
                for ( i = k + 1; i <= n - 1; i++) { //8
                       
                        
                         start_index7 = getIndex(i*(-1));

                       
                         start_index8 = getIndex(i-k);
                        int length1 = n - (i - k);
                        length1 = length1 - k;
                        //float start4 = clock();
                        for ( j = 0; j <= length1 - 1; j++) { //9
                                
                                 test =  C[start_index + i + j - k] + diagMatrixA[start_index7] * diagMatrixB[start_index8];
                                C[start_index + i + j - k] = test;
                                start_index7++;
                                start_index8++;

                        } //9

                        //float end4 = clock();
                        //total_time = total_time + (end4 - start4);

                        
                         start_index9 = getIndex(i-k);

                       
                        start_index10 = getIndex(i*(-1));
                      //  float start5 = clock();

                        for ( j = 0; j <= length1 - 1; j++) { //10

                                test = C[start_index + j] + diagMatrixA[start_index9 + k] * diagMatrixB[start_index10];

                                C[start_index + j] =test ;
                                start_index9++;
                                start_index10++;

                        } //10

                       // float end5 = clock();

                       // total_time = total_time + (end5 - start5);

                } //8

                for ( i = 0; i <= k; i++) { //11
                       
                        int temp6 = i;
                         start_index11 = 0;
                        if (temp6 == 0) {
                                start_index11 = temp6 * n - temp6 * (temp6 - 1) / 2;
                        } else {
                                start_index11 = n * (n + 1) / 2 + (temp6 - 1) * n - temp6 * (temp6 - 1) / 2;
                        }

                         start_index12 = 0;
                        
                        int temp7 = i - k;
                        int temp8 = abs(i - k);
                        if (temp7 < 0) {
                                start_index12 = n * (n + 1) / 2 + (temp8 - 1) * n - temp8 * (temp8 - 1) / 2;

                        } else {
                                start_index12 = temp7 * n - temp7 * (temp7 - 1) / 2;
                        }

                        
                        int length3 = n - i;
                        length3 = length3 - (k - i);

                       // float start6 = clock();

                        for ( j = 0; j <= length3 - 1; j++) { //12
 
                                test =  C[start_index + j] + diagMatrixA[start_index11 + (k - i)] * diagMatrixB[start_index12];
                                C[start_index + j] = test;
                                start_index11++;
                                start_index12++;

                        } //12
                      //  float end6 = clock();
                       // total_time = total_time + (end6 - start6);

                } //11

        } //7
  } //parallel

  return C;

}


 int Diagonal_Multiplication::getIndex(int k)
 {

 	 int abs_k;
	     int start_index;
		    if(k>=0)

		    {
		        start_index=k*dim-k*(k - 1) / 2;
	 
		    }
	
		    else
		    {
		       abs_k=abs(k);
		       start_index=dim * (dim + 1) / 2 + (abs_k - 1) * dim - abs_k * (abs_k - 1) / 2;

		    }

	     return start_index;
 }













