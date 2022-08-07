#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <stdio.h>

// x is the starting x; div is the 1d length 
void save_partial_matrix(int world_size, int step, int N, int x, int y, int rank, int div, double** C){
    char filename[50];
    sprintf(filename, "matrix_%d_%d.txt", rank, step);
    FILE* fh = fopen(filename, "w");
    fprintf(fh, "%d %d %d %d\n",N , x, y, div);
    for(int i = 1; i <= div; i++){
        for (int j = 1; j <= div; j++){
            fprintf(fh, "%e ", C[i][j]);
        }
        fprintf(fh, "\n");
    }
    fclose(fh);
}




int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //local index should be within [1, div] or [0, div + 1]
    //global index x will be local index + div * (rank / k)
    //global index y will be local index + div * (rank % k)
    int S = 0;
    int N, NT;
    double L, T, u, v;
    double dx, dy, dt;
    int k = sqrt(world_size);
    if(rank == 0){//read input parameters
    /*
        • N = 400 (Matrix Dimension)
        • NT = 20000 (Number of timesteps)
        • L = 1.0 (Physical Cartesian Domain Length)
        • T = 1.0e6 (Total Physical Timespan)
        • u = 5.0e-7 (X velocity Scalar)
        • v = 2.85e-7 (Y velocity Scalar)
        • S = 500 (save matrix per S steps)
    */
        N = atoi(argv[1]);
        printf("N = %d\n", N);

        NT = atoi(argv[2]);
        printf("NT = %d\n", NT);

        L = atof(argv[3]);
        printf("L = %e\n", L);

        T = atof(argv[4]);
        printf("T = %e\n", T);

        u = atof(argv[5]);
        printf("u = %e\n", u);
        
        v = atof(argv[6]);
        printf("v = %e\n", v);

        if(argc == 8){
            S = atoi(argv[7]);
            printf("S = %d\n", S);
        }

        dx = L / N, dy = L / N, dt = T / NT;
        if (dt > dx / sqrt(2 * (u * u + v * v))) {
            printf("%s", "Courant stability condition not satisfied!\n");
            return 0;
        }

        int *send_int_buf = malloc(sizeof(int) * 3);
        for(int rank = 1; rank < world_size; rank++) {
            send_int_buf[0] = N;
            send_int_buf[1] = NT;
            send_int_buf[2] = S;
            MPI_Send(send_int_buf, 3, MPI_INT, rank, 0, MPI_COMM_WORLD);
        }
        
        double *send_double_buf = malloc(sizeof(double) * 7);
        for(int rank = 1; rank < world_size; rank++) {
            send_double_buf[0] = L;
            send_double_buf[1] = T;
            send_double_buf[2] = u;
            send_double_buf[3] = v;
            send_double_buf[4] = dx;
            send_double_buf[5] = dy;
            send_double_buf[6] = dt;
            MPI_Send(send_double_buf, 7, MPI_DOUBLE, rank, 1, MPI_COMM_WORLD);
        }

        // Node 0 also needs to do its own calculation.
    } else {
        int *recv_int_buf = (int*) malloc(sizeof(int) * 3);
        MPI_Recv(recv_int_buf, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        N = recv_int_buf[0];
        NT = recv_int_buf[1];
        S = recv_int_buf[2];

        double *recv_double_buf = (double*) malloc(sizeof(double) * 7);
        MPI_Recv(recv_double_buf, 7, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        L = recv_double_buf[0];
        T = recv_double_buf[1];
        u = recv_double_buf[2];
        v = recv_double_buf[3];
        dx = recv_double_buf[4];
        dy = recv_double_buf[5];
        dt = recv_double_buf[6];
    }

    // for every node:
    int div = N/k;
    int length = div + 2;
    double **C = malloc(sizeof(double*) * length);
    double **Cn = malloc(sizeof(double*) * length);

    for (int i = 0; i < length; i++) {
        C[i] = malloc(sizeof(double) * length);
        Cn[i] = malloc(sizeof(double) * length);
    }

    for (int i = 0; i < length; i++){
        C[i][0] = 0;//left
        C[i][div + 1] = 0;//right
        C[0][i] = 0;//up
        C[div + 1][i] = 0;//down
    }
    
    // init starts
    double denom = 2 * (L / 4) * (L / 4);
    double x0 = L / 2, y0 = L / 2;//x0 and y0 are global idx

    // #pragma omp parallel for default(none) shared(C, denom, x0, y0, dx, dy, N)

    //i and j are local idx
    // x and y are global
    int offset_x = div * (rank / k);
    int offset_y = div * (rank % k);
    for (int i = 1; i <= div; i++) {
        for (int j = 1; j <= div; j++) {
            double x = (i - 1 + offset_x) * dx;
            double y = (j - 1 + offset_y) * dy;
            C[i][j] = exp(-(pow(x - x0, 2) + pow(y - y0, 2)) / denom);
        }
    }


    // save_partial_matrix(world_size, 0, N, offset_x, offset_y, rank, div,C);

    // init ends
    if(rank == 0)
        printf("thread: %d\n", omp_get_num_procs());
    
    double* left_recv_buf = (double*) malloc(sizeof(double) * (div));
    double* left_send_buf = (double*) malloc(sizeof(double) * (div));
    double* right_recv_buf = (double*) malloc(sizeof(double) * (div));
    double* right_send_buf = (double*) malloc(sizeof(double) * (div));
    
    int up = rank < k ? rank - k + world_size : rank - k;
    int down = rank >= world_size - k ? rank + k - world_size : rank + k;
    int left = rank % k == 0 ? rank - 1 + k : rank - 1;
    int right = (rank + 1) % k == 0 ? rank + 1 - k : rank + 1;

    save_partial_matrix(world_size, 0, N, offset_x, offset_y, rank, div, C);

    double time_messaging = 0;
    double time_start = omp_get_wtime();
    for (int n = 1; n <= NT; n++) {
        if (S && n % S == 0) {
            save_partial_matrix(world_size, n, N, offset_x, offset_y, rank, div, C);
        }
        //send lower boarder
        double msg_start = omp_get_wtime();
        MPI_Send(&C[div][1], div, MPI_DOUBLE, down, 2,MPI_COMM_WORLD);
        MPI_Recv(&C[0][1], div, MPI_DOUBLE, up, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        //send upper boarder. 
        MPI_Send(&C[1][1], div, MPI_DOUBLE, up,3, MPI_COMM_WORLD);
        MPI_Recv(&C[div + 1][1], div, MPI_DOUBLE, down, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //send right boarder
        for (int i = 1; i <= div; i++) {
            right_send_buf[i-1] = C[i][div];
        }
        MPI_Send(right_send_buf, div, MPI_DOUBLE, right,4, MPI_COMM_WORLD);

        MPI_Recv(left_recv_buf, div, MPI_DOUBLE, left, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 1; i <= div; i++) {
            C[i][0] = left_recv_buf[i-1];
        }

        //send left boarder
        for (int i = 1; i <= div; i++) {
            left_send_buf[i-1] = C[i][1];
        }
        MPI_Send(left_send_buf, div, MPI_DOUBLE, left, 5, MPI_COMM_WORLD);

        MPI_Recv(right_recv_buf, length - 2, MPI_DOUBLE, right, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for (int i = 1; i < length-1; i++) {
            C[i][div+1] = right_recv_buf[i-1];
        }

        double msg_end = omp_get_wtime();
        time_messaging += (msg_end - msg_start);

        #pragma omp parallel for default(none) shared(C, Cn, div, dx, dt, u, v)
        for (int i = 1; i <= div; i++) {
            for (int j = 1; j <= div; j++) {
                double u_term = u >= 0 ? u * (C[i][j] - C[i-1][j]) : u * (C[i+1][j] - C[i][j]);
                double v_term = v >= 0 ? v * (C[i][j] - C[i][j-1]) : v * (C[i][j+1] - C[i][j]);
                Cn[i][j] = C[i][j] - dt / dx * (u_term + v_term);
            }
        }

        //swap two array
        double **tmp = Cn;
        Cn = C;
        C = tmp;
            
    }
    double time_end = omp_get_wtime();
    double time_elapsed = time_end - time_start;
    if(rank == 0){
        printf("Time taken in message passing: %.2fs\n", time_messaging);
        printf("Time taken: %.2fs\n", time_elapsed);
        printf("Grind Rate: %d step/s\n", (int) (NT / time_elapsed));
    }
    MPI_Finalize();
}