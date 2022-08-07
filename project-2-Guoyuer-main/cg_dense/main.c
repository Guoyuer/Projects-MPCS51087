#include "cg_dense.h"

int main(int argc, char *argv[]) {
    // n is the physical domain size;
    int n = 100;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    printf("Solving Poisson Equation on %d x %d domain...\n", n, n);
    // Dimension of operator matrix and vectors is n^2
    int N = n * n;
    // Allocate full A matrix and vectors
    printf("Dense Memory  = %.2lf MB\n", (N * N + 5 * N) * sizeof(double) / 1024.0 / 1024.0);
//    return 0;
    double **A = matrix(N);
    double *x = (double *) calloc(N, sizeof(double));
    double *b = (double *) calloc(N, sizeof(double));

    // Compute elements of 'A' matrix (Poisson Operator)
    fill_A(A, N);

    // Compute elements of boundary condition vector 'b'
    fill_b(b, N);

    // Run Dense CG Solve
    double start = get_time();
    cg_dense(A, x, b, N);
    double stop = get_time();
    printf("Dense Runtime = %.2lf seconds\n", stop - start);

    // Save Solution Vector to File
    save_vector(x, N, "dense.out");

    // Free A matrix
    matrix_free(A);

    // Free vectors
    free(x);
    free(b);
}
