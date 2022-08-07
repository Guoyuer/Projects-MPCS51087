#include "cg_sparse_mpi.h"

// Parallel Distributed MPI Version of
// Specific "On the fly" Matrix Vector Product for v = A * w
// where 'A' is the 2D Poisson operator matrix
// that is applied without explicit storage
// v, w = 1D vectors
// N = dimension
void parallel_matvec_OTF( double * v, double * w, long N, int mype, int nprocs, MPI_Comm comm1d){
    long Nm = N / nprocs;
    long n = sqrt(N);
    int left, right; //adjacent ranks
    MPI_Cart_shift(comm1d, 0, 1, &left, &right);
    MPI_Request request[4];
    MPI_Status status[4];
    //if src or dst is MPI_PROC_NULL, recv or send will become dummy
    double *right_ghost = calloc(n, sizeof(double));
    double *left_ghost = calloc(n, sizeof(double));
    //send my left to left rank
    MPI_Isend(&w[0], n, MPI_DOUBLE, left, 99, comm1d, &request[0]);
    //recv my right from right rank
    MPI_Irecv(right_ghost, n, MPI_DOUBLE, right, 99, comm1d, &request[1]);

    //send my right to right rank
    MPI_Isend(&w[Nm - n], n, MPI_DOUBLE, right, 99, comm1d, &request[2]);
    MPI_Irecv(left_ghost, n, MPI_DOUBLE, left, 99, comm1d, &request[3]);
    MPI_Waitall(4, request, status);

    for(long i = 0; i < Nm; i++){
        long global_idx = i + Nm * mype;
        double left_val = 0, far_left_val = 0;
        double right_val = 0, far_right_val = 0;

        if(global_idx % n != 0){
            left_val = i == 0 ? left_ghost[n-1] : w[i-1];
        }

        if(global_idx % n != n-1){
            right_val = i == Nm-1 ? right_ghost[0] : w[i+1];
        }

        far_left_val = i - n < 0 ? left_ghost[i]: w[i-n];
        far_right_val = i + n >=Nm ? right_ghost[n - Nm + i] : w[i+n];

        v[i] = -far_left_val - left_val + 4 * w[i] - right_val - far_right_val;
    }

    

}



// Serial Conjugate Gradient Solver Function for Ax = b
// No need to store A. generate on the fly
// x = 1D solution vector
// b = 1D vector
// N = dimension
void parallel_cg_sparse_poisson(double * x, double * b, long N, int mype, int nprocs, MPI_Comm comm1d, int save_output){
    //x,b,r,p,z are all local
    long Nm = N / nprocs;
    double *r = (double *) malloc(Nm * sizeof(double)); //r is local

    // r = -A*x + b
    parallel_matvec_OTF(r, x, N, mype, nprocs, comm1d);

    parallel_axpy(-1.0, r, 1.0, b, N, mype, nprocs);

    //p = r;
    double *p = (double *) malloc(Nm * sizeof(double));
    memcpy(p, r, Nm * sizeof(double));

    //rsold = r' * r;
    double rsold = parallel_dotp(r, r, N, mype, nprocs, comm1d);
    
    double *z = (double *) malloc(Nm * sizeof(double));

    long iter = 0;
    for (iter = 0; iter < N; iter++) {
        //z = A * p;
        parallel_matvec_OTF(z, p, N, mype, nprocs, comm1d);

        //alpha = rsold / (p' * z);
        double alpha = rsold / parallel_dotp(p, z, N, mype, nprocs, comm1d);

        //x = x + alpha * p;
        parallel_axpy(1.0, x, alpha, p, N, mype, nprocs);

        //r = r - alpha * z;
        parallel_axpy(1.0, r, -alpha, z, N, mype, nprocs);

        double rsnew = parallel_dotp(r, r, N, mype, nprocs, comm1d);

        if (sqrt(rsnew) < 1.0e-10){
            if(save_output){
                if(mype == 0){
                    double *global_x = (double *) calloc(N, sizeof(double));
                    MPI_Gather(x, Nm, MPI_DOUBLE, global_x, Nm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                    save_vector(global_x, N, "global_x_mpi.out");
                }else{
                    MPI_Gather(x, Nm, MPI_DOUBLE, NULL    , 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                }
            }
            break;
        }

        //p = (rsnew / rsold) * p + r;
        parallel_axpy(rsnew / rsold, p, 1.0, r, N, mype, nprocs);

        rsold = rsnew;
    }
    if(mype == 0){
        printf("CG converged in %ld iterations.\n", iter);
    }

    free(r);
    free(p);
    free(z);
}


// Parallel Distributed MPI Version of
// Dot product of c = a * b
// c = result scalar that's returned
// a, b = 1D Vectors
// N = dimension
double parallel_dotp( double * a, double * b, long N, int mype, int nprocs, MPI_Comm comm1d)
{

    double local_sum = 0.0;
    for (long i = 0; i < N / nprocs; i++){
        local_sum += a[i] * b[i];
    }

    double global_sum = 0.0;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm1d);
    //every rank gets correct result.
    return global_sum;
}


// Parallel Distributed MPI Version of
// Scale and add of two vectors (axpy)
// Solves w = alpha * w + beta * v (overwrites what's in w)
// alpha, beta = scalars
// w, v = 1D vectors
// N = dimension
// every time the function gets called, w and v are local
void parallel_axpy( double alpha, double * w, double beta, double * v, long N, int mype, int nprocs)
{
    for (long i = 0; i < N / nprocs; i++)
        w[i] = alpha * w[i] + beta * v[i];
}


// Parallel Distributed MPI Version of
// Fills a 1-D RHS vector specifying boundary conditions
// by calling the get_b method
// b = 1-D RHS vector
// N = dimension (length of b)
void parallel_fill_b(double * b, long N, int mype, int nprocs){
    long n = sqrt(N);
    for(long idx = 0; idx < N / nprocs; idx++){
        long i = idx / n + n / nprocs * mype;
        long j = idx % n;
        b[idx] = find_b(i, j, n);
    }
}

// Sets a circular source term for the right hand side 'b' vector
// Takes as input 2D spatial domain indices
// i, j = indices
// n = physical dimension
double find_b(long i, long j, long n) {
    double delta = 1.0 / (double) n;

    double x = -.5 + delta + delta * j;
    double y = -.5 + delta + delta * i;

    // Check if within a circle
    double radius = 0.1;
    if (x * x + y * y < radius * radius)
        return delta * delta / 1.075271758e-02;
    else
        return 0.0;
}