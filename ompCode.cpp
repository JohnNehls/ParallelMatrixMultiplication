#include<stdio.h>
#include<cstdlib>
#include<time.h>
#include<iostream>
#include<matrix.h>
#include<omp.h>


/* Function Prototypes */
void FillMatricesRandomly(Matrix<double> &mat_a, Matrix<double> &mat_b ); /* fill the matricies RANDOMLY */
void PrintMatrices(Matrix<double> &mat_a, Matrix<double> &mat_b );	    /* print the input matricies and resultant matrix to the screen */

/* Instantiate Global Variables */
double start_time, end_time; 	/* used to time the multiplication */
int row_start, row_end; 	/* which rows of mat_a that are calculated by the slave process */
int randomHigh = 100; 		/* the upper bound of the random numbers that fill mat_a and mat_b*/
int randomLow = 0; 		/* the lower bound of the random numbers that fill mat_a and mat_b*/


int main(int argc, char *argv[]) {
  if ( argv[1]== NULL || argv[2] == NULL){
   std::cout << "ERROR: The program must be executed in the following way  \n\n  \t \"./omp.exe NumberOfThreads N \"  \n\n where NuberOfThreads and N are integers. \n \n " << std::endl;
   return 1;
  }

  int numThreads = atoi(argv[1]);
  std::cout << "The number of OpenMP threads: " << numThreads << std::endl;
  omp_set_dynamic(0);
  omp_set_num_threads(numThreads);

  int N = atoi(argv[2]); 	// for simplicity, all matricies will be NxN
  int NUM_ROWS_A = N;
  int NUM_COLUMNS_A = N;
  int NUM_ROWS_B = N;
  int NUM_COLUMNS_B = N;
  

  
  struct timespec start, end;

  /* Before anything else, check if matricies can be multiplied */
  if ( NUM_COLUMNS_A != NUM_ROWS_B){
    std::cout << "Program ended without a result: the matricies must be of the form mxn * nxp, where the inner dimensions agree." << std::endl;
      return 1;
	}
  Matrix<double> mat_a = Matrix<double>(NUM_ROWS_A, NUM_COLUMNS_A);
  Matrix<double> mat_b = Matrix<double>(NUM_ROWS_B, NUM_COLUMNS_B);
  Matrix<double> mat_result = Matrix<double>(NUM_ROWS_A, NUM_COLUMNS_B);
  
  FillMatricesRandomly(mat_a, mat_b); 	/* Excluded from timing!!!  */

 clock_gettime(CLOCK_MONOTONIC, &start);
  
  // Matrix Multiplication
#pragma omp parallel for 
    for (int i = 0; i < NUM_ROWS_A; i++) {
      for (int j = 0; j < NUM_COLUMNS_B; j++) {
  	for (int k = 0; k < NUM_ROWS_B; k++) {
  	  mat_result(i,j) += ( mat_a(i,k) * mat_b(k,j) );
  	}
      }
    }

    clock_gettime(CLOCK_MONOTONIC, &end);
    double diffTime = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)*1e-9;
    std::cout << "\nRunning Time = " << diffTime << std::endl;
    // PrintMatrices(mat_a, mat_b);
  return 0;
}
void FillMatricesRandomly(Matrix<double> &mat_a, Matrix<double> &mat_b){
  srand( time( NULL ));		/* initialize the random number generator with the current time */

  for (int i = 0; i < mat_a.rows(); i++) {
    for (int j = 0; j < mat_a.cols(); j++) {
      mat_a(i,j) = rand() % (randomHigh - randomLow) + randomLow;      
    }
  }
  for (int i = 0; i < mat_a.rows(); i++) {
    for (int j = 0; j < mat_a.cols(); j++) {
      mat_b(i,j) = rand() % (randomHigh - randomLow) + randomLow;      
    }
  }
}

void PrintMatrices(Matrix<double> &mat_a, Matrix<double> &mat_b){

  for (int i = 0; i < mat_a.rows(); i++) {
    printf("\n");
    for (int j = 0; j < mat_b.cols(); j++)
      std::cout << mat_a(i,j) << " ";
  }
  printf("\n");

  for (int i = 0; i < mat_b.rows(); i++) {
    printf("\n");
    for (int j = 0; j < mat_b.cols(); j++)
      std::cout << mat_b(i,j) << " ";
  }
  printf("\n"); 
    }







