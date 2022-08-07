hybrid: main.c
	mpicc -O2 -fopenmp -lm main.c -o hybrid_adv
first_order: first_order.c
	mpicc -O2 -fopenmp -lm first_order.c -o first_order
second_order: second_order.c
	mpicc -O2 -fopenmp -lm second_order.c -o second_order
test: test.cpp
	mpicxx test.cpp -o test