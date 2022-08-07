#include "cg_sparse_mpi.h"

#define DIMENSION 1

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    // n is the physical domain size;
    int n = 100;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    // Dimension of operator matrix and vectors is n^2
    int N = n * n;
    int Nm = N / nprocs;
    if(mype == 0){
        printf("Solving Poisson Equation on %d x %d domain...\n", n, n);
        printf("Sparse Memory  = %.2lf MB\n", (5 * N) * sizeof(double) / 1024.0 / 1024.0);
    }
    double *x = (double *) calloc(Nm, sizeof(double));
    double *b = (double *) calloc(Nm, sizeof(double));

    parallel_fill_b(b, N, mype, nprocs);

    // //global_b is for debug purposes
    // double *global_b;
    // if(mype == 0){
    //     global_b = (double *) calloc(N, sizeof(double));
    //     MPI_Gather(b, Nm, MPI_DOUBLE, global_b, Nm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //     save_vector(global_b, N, "sparse_b.out");
    // }else{
    //     MPI_Gather(b, Nm, MPI_DOUBLE, NULL    , 0 , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // }



    // Run Dense CG Solve
    // // global_x is used to store concatenated local_x
    
    int dims[DIMENSION], periodic[DIMENSION];
	MPI_Comm comm1d;
	dims[0] = nprocs;    // Number of MPI ranks in each dimension
	periodic[0] = 0;             // Turn on/off periodic boundary conditions for each dimension

	// Create Cartesian Communicator
	MPI_Cart_create( MPI_COMM_WORLD, // Starting communicator we are going to draw from
			DIMENSION,      // MPI grid n-dimensionality
			dims,           // Array holding number of MPI ranks in each dimension
			periodic,       // Array indicating if we want to use periodic BC's or not
			1,              // Yes/no to reordering (allows MPI to re-organize for better perf)
			&comm1d );      // Pointer to our new Cartesian Communicator object

    int save_output = 0;
    double start = get_time();
    parallel_cg_sparse_poisson(x, b, N, mype, nprocs, comm1d, save_output);
    double stop = get_time();
    if(mype == 0){
        printf("nodes = %d\n", nprocs);
        printf("Runtime = %.2lf seconds\n", stop - start);
    }

    // Free vectors
    free(x);
    free(b);
    MPI_Finalize();
    return 0;
}