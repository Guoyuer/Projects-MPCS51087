#ifndef ch_sparse_mpi_h
#define ch_sparse_mpi_h


#endif /* ch_sparse_mpi_h */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>


#ifdef MPI
#endif

void parallel_cg_sparse_poisson(double * x, double * b, long N, int mype, int nprocs,MPI_Comm comm1d, int save_output);
void parallel_matvec_OTF( double * v, double * w, long N, int mype, int nprocs, MPI_Comm comm1d);
void parallel_axpy( double alpha, double * w, double beta, double * v, long N, int mype, int nprocs);
double parallel_dotp( double * a, double * b, long N, int mype, int nprocs, MPI_Comm comm1d);
double find_b(long i, long j, long n);
void parallel_fill_b(double * b, long N, int mype, int nprocs);


// utils.c
void print_matrix(double ** A, long N );
void print_vector(double * x, long N );
void save_vector(double * x, long N, char * fname );
double ** matrix( long N );
void matrix_free( double ** M);
double get_time(void);
void cli_error(void);

