#include "cg_sparse.h"

// Specific "On the fly" Matrix Vector Product for v = A * w
// where 'A' is the 2D Poisson operator matrix
// that is applied without explicit storage
// v, w = 1D vectors
// N = dimension
void matvec_fly(double *v, double *w, long N) {
    // Set solution vector to 0
    long n = sqrt(N);
    memset(v, 0, N * sizeof(double));
//  v[i] = −w[i − n] − w[i − 1] + 4w[i] − w[i + 1] − w[i + n]
//  v[i] = -t1       - t2       + t3    - t4       - t5
//          far right    right      center   left       far left?
    for (long i = 0; i < N; i++) {
        long x = i % n;
        double t1, t2, t3, t4, t5;
        t1 = i - n >= 0 ? w[i - n] : 0;
        t2 = i - 1 >= 0 && x != 0 ? w[i - 1] : 0;
        t3 = 4 * w[i];
        t4 = i + 1 < N && x != n - 1 ? w[i + 1] : 0;
        t5 = i + n < N ? w[i + n] : 0;
        v[i] = -t1 - t2 + t3 - t4 - t5;
    }
}


// Serial Conjugate Gradient Solver Function for Ax = b
// No need to store A. generate on the fly
// x = 1D solution vector
// b = 1D vector
// N = dimension
void cg_sparse(double *x, double *b, long N) {
    // r = -A*x + b
    double *r = (double *) malloc(N * sizeof(double));
    matvec_fly(r, x, N);
    axpy(-1.0, r, 1.0, b, N);

    //p = r;
    double *p = (double *) malloc(N * sizeof(double));
    memcpy(p, r, N * sizeof(double));

    //rsold = r' * r;
    double rsold = dotp(r, r, N);

    // z
    double *z = (double *) malloc(N * sizeof(double));

    long iter = 0;
    for (iter = 0; iter < N; iter++) {
        //z = A * p;
        matvec_fly(z, p, N);

        //alpha = rsold / (p' * z);
        double alpha = rsold / dotp(p, z, N);

        //x = x + alpha * p;
        axpy(1.0, x, alpha, p, N);

        //r = r - alpha * z;
        axpy(1.0, r, -alpha, z, N);

        double rsnew = dotp(r, r, N);

        if (sqrt(rsnew) < 1.0e-10)
            break;

        //p = (rsnew / rsold) * p + r;
        axpy(rsnew / rsold, p, 1.0, r, N);

        rsold = rsnew;
    }
    printf("CG converged in %ld iterations.\n", iter);

    free(r);
    free(p);
    free(z);
}


// Dot product of c = a * b
// c = result scalar that's returned
// a, b = 1D Vectors
// N = dimension
double dotp(double *a, double *b, long N) {
    double c = 0.0;
    for (long i = 0; i < N; i++)
        c += a[i] * b[i];

    return c;
}

// Scale and add of two vectors (axpy)
// Solves w = alpha * w + beta * v (overwrites what's in w)
// alpha, beta = scalars
// w, v = 1D vectors
// N = dimension
void axpy(double alpha, double *w, double beta, double *v, long N) {
    for (long i = 0; i < N; i++)
        w[i] = alpha * w[i] + beta * v[i];
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

// Fills a 1-D RHS vector specifying boundary conditions
// by calling the get_b method
// b = 1-D RHS vector
// N = dimension (length of b)
void fill_b(double *b, long N) {
    long n = sqrt(N);
    for (long i = 0; i < n; i++)
        for (long j = 0; j < n; j++) {
            long idx = i * n + j;
            b[idx] = find_b(i, j, n);
        }
}