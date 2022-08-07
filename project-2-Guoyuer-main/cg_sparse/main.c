#include "cg_sparse.h"


int main(int argc, char *argv[]) {
    // n is the physical domain size;
    int n = 100;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    printf("Solving Poisson Equation on %d x %d domain...\n", n, n);
    // Dimension of operator matrix and vectors is n^2
    int N = n * n;
    printf("Sparse Memory  = %.2lf MB\n", (5 * N) * sizeof(double) / 1024.0 / 1024.0);
//    return 0;
    double *x = (double *) calloc(N, sizeof(double));
    double *b = (double *) calloc(N, sizeof(double));

    // Compute elements of boundary condition vector 'b'
    fill_b(b, N);

    save_vector(b, N, "sparse_b_ref.out");

    // Run Dense CG Solve
    double start = get_time();
    cg_sparse(x, b, N);
    double stop = get_time();
    printf("Sparse Runtime = %.2lf seconds\n", stop - start);

    // Save Solution Vector to File
    save_vector(x, N, "sparse.out");

    // Free vectors
    free(x);
    free(b);
}