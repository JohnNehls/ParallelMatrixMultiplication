#include<cstdlib>
#include<time.h>
#include<iostream>
#include<matrix.h>
#include<omp.h>


/* Function Prototypes */
void FillMatricesRandomly(Matrix<double> &A, Matrix<double> &B );
void PrintMatrices( Matrix<double> &A, Matrix<double> &B, Matrix<double> &C );

/* Global variables that could be inputs #TODO */
int randomHigh = 100; 		/* the upper bound of the random numbers that fill A and B*/
int randomLow = 0; 		/* the lower bound of the random numbers that fill A and B*/

int main(int argc, char *argv[]) {
    std::cout << "Starting an OpenMP parallel matrix multiplication. \n  " << std::endl;  
  // Read in the two inputs: NumberOfOMPthreads and N
  if ( argv[1]== NULL || argv[2] == NULL){ // check if inputs were supplied
    std::cout << "ERROR: The program must be executed in the following way  \n\n  \t \"./omp.exe NumberOfThreads N \"  \n\n where NuberOfThreads and N are integers. \n \n " << std::endl;
    return 1;
  }
  int numThreads = atoi(argv[1]); 
  std::cout << "The number of OpenMP threads: " << numThreads << std::endl;
  omp_set_dynamic(0);		// do not allow the number of threads to be set internally
  omp_set_num_threads(numThreads); // set the number of threads 

  int N = atoi(argv[2]);  // The dimensions of the matrices MUST be specified at runtime.
  std::cout << "The matrices are: " << N<<"x"<<N<< std::endl;
  // for simplicity, all matricies will be NxN
  int numberOfRowsA = N;  int numberOfColsA = N;  int numberOfRowsB = N;  int numberOfColsB = N;  
  //Declare matrices: Matix class is a 2D vetor
  Matrix<double> A = Matrix<double>(numberOfRowsA, numberOfColsA); 
  Matrix<double> B = Matrix<double>(numberOfRowsB, numberOfColsB); 
  Matrix<double> C = Matrix<double>(numberOfRowsA, numberOfColsB); 

  FillMatricesRandomly(A, B); 	/* Excluded from timing!!!  */

  struct timespec start, end;
  clock_gettime(CLOCK_MONOTONIC, &start); // start timing

  // Used to calculate the longest a process spends calculating its part of the workload.
  // double sumLocalTime = 0;
  double sum = 0;
  double val = 0;
  // Matrix Multiplication
#pragma omp parallel for private(val) reduction(+:sum)
  for (int i = 0; i < A.rows(); i++) {//iterate through rows of A (parallelized loop)
    val = omp_get_wtime();
    for (int j = 0; j < B.cols(); j++) {//iterate through columns of B
      for (int k = 0; k < B.rows(); k++) {//iterate through rows of B
	C(i,j) += (A(i,k) * B(k,j));
      }
    }
    sum +=  omp_get_wtime() - val;
  }

  clock_gettime(CLOCK_MONOTONIC, &end); // end timing
  double totalMatrixCalculationTime = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)*1e-9;
  std::cout << "Total multplication time = " << totalMatrixCalculationTime << std::endl;
  std::cout << "average multplication time = " << sum/numThreads << std::endl;
  std::cout << "Approximate Communication time = " << totalMatrixCalculationTime - sum/numThreads << std::endl;  
//   PrintMatrices(A, B, C); 
  return 0;
}

void FillMatricesRandomly(Matrix<double> &A, Matrix<double> &B){
/* initialize the random number generator with the current time */  
  srand( time( NULL ));		
  for (int i = 0; i < A.rows(); i++) {
    for (int j = 0; j < A.cols(); j++) {
      A(i,j) = rand() % (randomHigh - randomLow) + randomLow;      
    }
  }
  for (int i = 0; i < B.rows(); i++) {
    for (int j = 0; j < B.cols(); j++) {
      B(i,j) = rand() % (randomHigh - randomLow) + randomLow;      
    }
  }
}

void PrintMatrices(Matrix<double> &A, Matrix<double> &B, Matrix<double> &C){
  for (int i = 0; i < A.rows(); i++) {
    std::cout <<"\n"<<std::endl;
    for (int j = 0; j < A.cols(); j++)
      std::cout << A(i,j) << " ";
  }
    std::cout <<"\n\n"<<std::endl;  
  for (int i = 0; i < B.rows(); i++) {
    std::cout <<"\n"<<std::endl;
    for (int j = 0; j < B.cols(); j++)
      std::cout << B(i,j) << " ";
  }
  std::cout <<"\n\n"<<std::endl;  
  for (int i = 0; i < C.rows(); i++) {
    std::cout <<"\n"<<std::endl;  
    for (int j = 0; j < C.cols(); j++)
      std::cout << C(i,j) << " ";
  }
  std::cout <<"\n\n"<<std::endl;    

}



 
