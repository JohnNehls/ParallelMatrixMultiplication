#include<stdio.h>
#include<cstdlib>
#include<ctime>
#include<iostream>
#include<mpi.h>
#include<matrix.h>

// Function Prototypes
void FillMatricesRandomly(Matrix<double> &mat_a, Matrix<double> &mat_b ); 
void printArrays( Matrix<double> &mat_a, Matrix<double> &mat_b, Matrix<double> &mat_result );


/* MPI Send and Recieve Tags */
#define ROW_START_TAG 0 //tag for messages sent from master to slaves
#define ROW_END_TAG 1 //tag for messages sent from master to slaves
#define MAT_A_ROWS_TAG 2 //tag for messages sent from master to slaves
#define MAT_RESULT_ROWS_TAG 3 //tag for messages sent from master to slaves

// Instantiate Global Variables
int rank;             // mpi: process id number
int nProcesses;      // mpi: number of total processess 
int row_start, row_end;       // which rows of mat_a that are calculated by the slave process
int granularity; 		// granularity of parallelization (# of rows per processor) 
MPI_Status status; // store status of a MPI_Recv
MPI_Request request; //capture request of a MPI_Isend

// Used to time the multiplication 
double start_time, end_time;

/* Defines that COULD be inputs oneday #TODO */
#define NUM_ROWS_A 240
#define NUM_COLUMNS_A 240
#define NUM_ROWS_B 240
#define NUM_COLUMNS_B 240
int randomHigh = 3; 		// the upper bound of the random numbers that fill mat_a and mat_b
int randomLow = 0; 		// the lower bound of the random numbers that fill mat_a and mat_b


int main(int argc, char *argv[]) {
  

  Matrix<double> mat_a = Matrix<double>(NUM_ROWS_A, NUM_COLUMNS_A);
  Matrix<double> mat_b = Matrix<double>(NUM_ROWS_B, NUM_COLUMNS_B);
  Matrix<double> mat_result = Matrix<double>(NUM_ROWS_A, NUM_COLUMNS_B);


  /* Before anything else, check if matricies can be multiplied */
  if ( NUM_COLUMNS_A != NUM_ROWS_B){
    std::cout << "Program ended without a result: the matricies must be of the form mxn * nxp, where the inner dimensions agree." << std::endl;
    return 1;
  }
    
  MPI_Init(&argc, &argv);		      /* initialize MPI */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);	      /* store the rank */
  MPI_Comm_size(MPI_COMM_WORLD, &nProcesses); /* store the number of processes */

  /* master initializes work*/
  if (rank == 0) {
    FillMatricesRandomly( mat_a, mat_b); 	/* Excluded from timing!!!  */    

    start_time = MPI_Wtime(); 	/* Begin Timing */    
    for (int i = 1; i < nProcesses; i++) { /* for each slave */

      // calculate granularity (-1 comes from excluding the master process)      
      granularity = (NUM_ROWS_A / (nProcesses - 1)); 
      row_start = (i - 1) * granularity;
      if (((i + 1) == nProcesses) && ((NUM_ROWS_A % (nProcesses - 1)) != 0)) {//if rows of [A] cannot be equally divided among slaves
	row_end = NUM_ROWS_A; //last slave gets all the remaining rows
      } else {
	row_end = row_start + granularity; //rows of [A] are equally divisable among slaves
      }
      //send the low bound, without blocking, to the intended slave 
      MPI_Isend(&row_start, 1, MPI_INT, i, ROW_END_TAG, MPI_COMM_WORLD, &request);
      //next send the upper bound without blocking, to the intended slave
      MPI_Isend(&row_end, 1, MPI_INT, i ,ROW_START_TAG, MPI_COMM_WORLD, &request);
      //finally send the allocated row granularity of [A] without blocking, to the intended slave
      MPI_Isend(&mat_a(row_start,0), (row_end - row_start) * NUM_COLUMNS_A, MPI_DOUBLE, i, MAT_A_ROWS_TAG, MPI_COMM_WORLD, &request);  
    }
  }
  //broadcast mat_b (MPI_Bcast: Broadcasts a message from the process with rank "root" to all other processes of the communicator)
  MPI_Bcast(&mat_b(0,0), NUM_ROWS_B*NUM_COLUMNS_B, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* work done by slaves (not rank = 0)*/
  if (rank > 0) {
    //receive low bound from the master
    MPI_Recv(&row_start, 1, MPI_INT, 0, ROW_END_TAG, MPI_COMM_WORLD, &status);
    //next receive upper bound from the master
    MPI_Recv(&row_end, 1, MPI_INT, 0, ROW_START_TAG, MPI_COMM_WORLD, &status);
    //finally receive row granularity of [A] to be processed from the master
    MPI_Recv(&mat_a(row_start,0), (row_end - row_start) * NUM_COLUMNS_A, MPI_DOUBLE, 0, MAT_A_ROWS_TAG, MPI_COMM_WORLD, &status);

    /* Matrix Multiplication  */
    for (int i = row_start; i < row_end; i++) {//iterate through a given set of rows of [A]
      for (int j = 0; j < NUM_COLUMNS_B; j++) {//iterate through columns of [B]
	for (int k = 0; k < NUM_ROWS_B; k++) {//iterate through rows of [B]
	  mat_result(i,j) += (mat_a(i,k) * mat_b(k,j));
	}
      }
    }
    
    //send back the low bound first without blocking, to the master
    MPI_Isend(&row_start, 1, MPI_INT, 0, ROW_END_TAG, MPI_COMM_WORLD, &request);
    //send the upper bound next without blocking, to the master
    MPI_Isend(&row_end, 1, MPI_INT, 0, ROW_START_TAG, MPI_COMM_WORLD, &request);
    //finally send the processed granularity of data without blocking, to the master
    MPI_Isend(&mat_result(row_start,0), (row_end - row_start) * NUM_COLUMNS_B, MPI_DOUBLE, 0, MAT_RESULT_ROWS_TAG, MPI_COMM_WORLD, &request);
    // MPI_Isend(&mat_result[row_start][0], (row_end - row_start) * NUM_COLUMNS_B, MPI_DOUBLE, 0, MAT_RESULT_ROWS_TAG, MPI_COMM_WORLD, &request);

  }
  /* master gathers processed work*/
  if (rank == 0) {
    for (int i = 1; i < nProcesses; i++) {// untill all slaves have handed back the processed data
      //receive low bound from a slave
      MPI_Recv(&row_start, 1, MPI_INT, i, ROW_END_TAG, MPI_COMM_WORLD, &status);
      //receive upper bound from a slave
      MPI_Recv(&row_end, 1, MPI_INT, i, ROW_START_TAG, MPI_COMM_WORLD, &status);
      // //receive processed data from a slave
      MPI_Recv(&mat_result(row_start,0), (row_end - row_start) * NUM_COLUMNS_B, MPI_DOUBLE, i, MAT_RESULT_ROWS_TAG, MPI_COMM_WORLD, &status);
      // MPI_Recv(&mat_result[row_start][0], (row_end - row_start) * NUM_COLUMNS_B, MPI_DOUBLE, i, MAT_RESULT_ROWS_TAG, MPI_COMM_WORLD, &status);

    }
    end_time = MPI_Wtime();
    printf("\nRunning Time = %f\n\n", end_time - start_time);
    // printArrays(mat_a,  mat_b,  mat_result);

  }
  MPI_Finalize(); //finalize MPI operations
  return 0;
}

void FillMatricesRandomly(Matrix<double> &mat_a, Matrix<double> &mat_b){
  srand( time( NULL ));		/* initialize the random number generator with the current time */

  for (int i = 0; i < NUM_ROWS_A; i++) {
    for (int j = 0; j < NUM_COLUMNS_A; j++) {
      mat_a(i,j) = rand() % (randomHigh - randomLow) + randomLow;      
    }
  }
  for (int i = 0; i < NUM_ROWS_B; i++) {
    for (int j = 0; j < NUM_COLUMNS_B; j++) {
      mat_b(i,j) = rand() % (randomHigh - randomLow) + randomLow;      
    }
  }
}

void printArrays(Matrix<double> &mat_a, Matrix<double> &mat_b, Matrix<double> &mat_result){

  for (int i = 0; i < NUM_ROWS_A; i++) {
    printf("\n");
    for (int j = 0; j < NUM_COLUMNS_A; j++)
      std::cout << mat_a(i,j) << " ";
  }
  printf("\n");

  for (int i = 0; i < NUM_ROWS_B; i++) {
    printf("\n");
    for (int j = 0; j < NUM_COLUMNS_B; j++)
      std::cout << mat_b(i,j) << " ";
  }
  printf("\n");
  for (int i = 0; i < NUM_ROWS_B; i++) {
    printf("\n");
    for (int j = 0; j < NUM_COLUMNS_B; j++)
      std::cout << mat_result(i,j) << " ";
  }
  printf("\n"); 

}

// void FillMatricesRandomly(){
//   srand( time( NULL ));		/* initialize the random number generator with the current time */

//   for (i = 0; i < NUM_ROWS_A; i++) {
//     for (j = 0; j < NUM_COLUMNS_A; j++) {
//       mat_a(i,j) = rand() % (randomHigh - randomLow) + randomLow;
//     }
//   }
//   for (i = 0; i < NUM_ROWS_B; i++) {
//     for (j = 0; j < NUM_COLUMNS_B; j++) {
//       mat_b(i,j) = rand() % (randomHigh - randomLow) + randomLow;
//       // mat_b[i][j] = rand() % (randomHigh - randomLow) + randomLow;
//     }
//   }
// }

// void FillMatrices(){
//   for (i = 0; i < NUM_ROWS_A; i++) {
//     for (j = 0; j < NUM_COLUMNS_A; j++) {
//       mat_a(i,j) = i + j;
//     }
//   }
//   for (i = 0; i < NUM_ROWS_B; i++) {
//     for (j = 0; j < NUM_COLUMNS_B; j++) {
//       mat_b(i,j) = i*j;
//     }
//   }
// }

// void printArrays(){
//   for (i = 0; i < NUM_ROWS_A; i++) {
//     printf("\n");
//     for (j = 0; j < NUM_COLUMNS_A; j++)
//       std::cout << mat_a(i,j) << " ";
//   }
//   printf("\n");
//   for (i = 0; i < NUM_ROWS_B; i++) {
//     printf("\n");
//     for (j = 0; j < NUM_COLUMNS_B; j++)
//       std::cout << mat_b(i,j) << " ";
//   }
//   printf("\n");
//   for (i = 0; i < NUM_ROWS_A; i++) {
//     printf("\n");
//     for (j = 0; j < NUM_COLUMNS_B; j++)
//       std::cout << mat_result[i][j] << " ";     
//   }
//   printf("\n");
// }
