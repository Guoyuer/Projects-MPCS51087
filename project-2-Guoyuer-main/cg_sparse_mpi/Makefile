mpi: main.c utils.c cg_sparse_mpi.c
	mpicc main.c utils.c cg_sparse_mpi.c -O2 -I . -lm -o main

test: test.c
	mpicc test.c -o test
