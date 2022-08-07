main_mpi: main_mpi.c
	mpicc main_mpi.c timer.c -o main_mpi -lm -O2 -fopenmp

mpi_test: mpi_test.c
	mpicc mpi_test.c -o mpi_test -lm -O2 -fopenmp


main: main.c
	gcc main.c timer.c -o main -lm -O2 -fopenmp