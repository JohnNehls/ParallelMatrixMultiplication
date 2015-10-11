#include<stdio.h>
#include<cstdlib>
#include<time.h>
#include<iostream>
#include<matrix.h>
#include<omp.h>

/* MPI Send and Recieve Tags */
#define ROW_START_TAG 0 //tag for messages sent from master to slaves
#define ROW_END_TAG 1 //tag for messages sent from master to slaves
#define MAT_A_ROWS_TAG 2 //tag for messages sent from master to slaves
#define MAT_RESULT_ROWS_TAG 3 //tag for messages sent from master to slaves

/* Function Prototypes */
void FillMatrices();	       /* fills the matricies to be multiplied */
void FillMatricesRandomly(Matrix<double> &mat_a, Matrix<double> &mat_b ); /* fill the matricies RANDOMLY */
void printArrays(Matrix<double> &mat_a, Matrix<double> &mat_b );	    /* print the input matricies and resultant matrix to the screen */

/* Instantiate Global Variables */
double start_time, end_time; 	/* used to time the multiplication */
int row_start, row_end; 	/* which rows of mat_a that are calculated by the slave process */
int randomHigh = 100; 		/* the upper bound of the random numbers that fill mat_a and mat_b*/
int randomLow = 0; 		/* the lower bound of the random numbers that fill mat_a and mat_b*/


int main(int argc, char *argv[]) {
  if ( argv[1]== NULL ){
   std::cout << "ERROR: The program must be executed in the following way  \n\n  \t \"./omp.exe NumberOfThreads N \"  \n\n where N is an integer. \n \n " << std::endl;
   return 1;
  }

  int N = atoi(argv[2]); 	// for simplicity, all matricies will be NxN
  int NUM_ROWS_A = N;
  int NUM_COLUMNS_A = N;
  int NUM_ROWS_B = N;
  int NUM_COLUMNS_B = N;
    
  struct timespec start, end;

  Matrix<double> mat_a = Matrix<double>(NUM_ROWS_A, NUM_COLUMNS_A);
  Matrix<double> mat_b = Matrix<double>(NUM_ROWS_B, NUM_COLUMNS_B);
  Matrix<double> mat_result = Matrix<double>(NUM_ROWS_A, NUM_COLUMNS_B);
  
  FillMatricesRandomly(mat_a, mat_b); 	/* Excluded from timing!!!  */

 clock_gettime(CLOCK_MONOTONIC, &start);
  
  // Matrix Multiplication
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
    // printArrays(mat_a, mat_b);
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

void printArrays(Matrix<double> &mat_a, Matrix<double> &mat_b){

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







