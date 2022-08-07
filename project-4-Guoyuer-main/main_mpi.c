#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include <time.h> //for random seed
#include <omp.h>
#include <stddef.h>
#include <mpi.h>
#include <string.h>

#define DIMENSION 1
#define SOFTENING 1e-3f
#define G 100
// ref: https://stackoverflow.com/questions/9864510/struct-serialization-in-c-and-transfer-over-mpi
typedef struct
{
    double m;          /* particle masses */
    double x, y, z;    /* particle positions */
    double vx, vy, vz; /* particle momenta */
} Particle;

double getRandIn(double min, double max)
{
    // 0 ~ 1 -> 0 ~ (max - min) -> min ~ max
    return (max - min) * ((rand() / (double)RAND_MAX)) + min;
}

// r: init for rank `r`
void init(Particle *p, int n, double L)
{
    L /= 2;
    for (int i = 0; i < n; i++)
    {
        p[i].m = 1.0f;
        p[i].x = getRandIn(-L, L);
        p[i].y = getRandIn(-L, L);
        p[i].z = getRandIn(-L, L); // on same plane: Z = 0
        p[i].vx = getRandIn(-10, 10);
        p[i].vy = getRandIn(-10, 10);
        p[i].vz = getRandIn(-10, 10);

        // p[i].m = 1.0f;

        // p[i].x = getRandIn(-5, 5);
        // p[i].y = getRandIn(-2, 2);
        // p[i].z = 0.0f; // on same plane: Z = 0
        // p[i].vx = 0.0f;
        // p[i].vy = 0.0f;
        // p[i].vz = 0.0f;
    }
}

/* calculate all interparticle forces and update instantaneous velocities */
void calc_f_and_update_v(Particle *local, Particle *remote, float dt, int nl)
{
#pragma omp parallel for
    for (int i = 0; i < nl; i++)
    {
        // printf("num of threads: %d\n", omp_get_num_threads());

        double Fx = 0.0f;
        double Fy = 0.0f;
        double Fz = 0.0f;

        for (int j = 0; j < nl; j++)
        {
            /* calculate net particle for on i'th particle */
            float dx = remote[j].x - local[i].x;
            float dy = remote[j].y - local[i].y;
            float dz = remote[j].z - local[i].z;
            float D_rij = sqrtf(dx * dx + dy * dy + dz * dz + SOFTENING); // distance of two bodies
            float D_rij3 = D_rij * D_rij * D_rij;
            float left_term = G * local[i].m * remote[j].m / D_rij3;

            Fx += dx * left_term;
            Fy += dy * left_term;
            Fz += dz * left_term;
        }
        /* update instantaneous velocity based on force and timestep */
        local[i].vx += dt * Fx;
        local[i].vy += dt * Fy;
        local[i].vz += dt * Fz;
    }
}

int main(int argc, char **argv)
{
    srand(time(NULL));
    int nParticles = atoi(argv[1]); /* number of particles */
    int nIters = atoi(argv[2]);     /* number of steps in simulation */
    int save = atoi(argv[3]);       /* if save == 1 then output to file*/
    double L = 20000;
    const float dt = 0.01f; /* time step   */
    MPI_Init(&argc, &argv);
    int M;
    MPI_Comm_size(MPI_COMM_WORLD, &M);
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    int nl = nParticles / M;

    MPI_Datatype MPI_BODY;
    MPI_Type_contiguous(7, MPI_DOUBLE, &MPI_BODY);
    MPI_Type_commit(&MPI_BODY);

    int dims[DIMENSION], periodic[DIMENSION];
    MPI_Comm comm1d;
    dims[0] = M;     // Number of MPI ranks in each dimension
    periodic[0] = 1; // Turn on/off periodic boundary conditions for each dimension
    // Create Cartesian Communicator
    MPI_Cart_create(MPI_COMM_WORLD, // Starting communicator we are going to draw from
                    DIMENSION,      // MPI grid n-dimensionality
                    dims,           // Array holding number of MPI ranks in each dimension
                    periodic,       // Array indicating if we want to use periodic BC's or not
                    1,              // Yes/no to reordering (allows MPI to re-organize for better perf)
                    &comm1d);       // Pointer to our new Cartesian Communicator object

    // local bodies
    Particle *B_local = malloc(nl * sizeof(Particle));
    // remote bodies
    Particle *B_remote = malloc(nl * sizeof(Particle));
    Particle *B_tmp_buf = malloc(nl * sizeof(Particle));

    // printf("world size: %d; mype: %d\n", M, mype);
    if (mype == 0)
    {
        // printf("nl: %d\n", nl);
        init(B_local, nl, L);
        for (int r = 1; r < M; r++) // only send to others.
        {
            init(B_tmp_buf, nl, L);
            MPI_Send(B_tmp_buf, nl, MPI_BODY, r, 99, comm1d);
        }
    }
    else
    {

        MPI_Recv(B_local, nl, MPI_BODY, 0, 99, comm1d, MPI_STATUS_IGNORE);
    }

    if (save)
    {
        if (mype == 0)
        {
            FILE *metadata = fopen("production_meta.dat", "w");
            fprintf(metadata, "%d %d\n", nParticles, nIters);
            fclose(metadata);
        }
    }
    MPI_Barrier(comm1d);
    /* ------------------------------*/
    /*     MAIN LOOP                 */
    /* ------------------------------*/

    if (mype == 0)
    {
        StartTimer();
    }
    MPI_Request request[2];
    MPI_Status status[2];
    for (int iter = 1; iter <= nIters; iter++)
    {
        if (mype == 0)
        {
            printf("iter: %d\n", iter);
        }
        memcpy(B_remote, B_local, nl * sizeof(Particle));
        FILE *step_file;
        if (save && mype == 0)
        {
            char filename[30];
            sprintf(filename, "production_step_%d.dat", iter);
            step_file = fopen(filename, "w");
        }
        for (int k = 0; k < M; k++)
        {
            calc_f_and_update_v(B_local, B_remote, dt, nl); /* compute interparticle forces and update vel */
            int left, right;                                // adjacent ranks
            MPI_Cart_shift(comm1d, 0, 1, &left, &right);
            MPI_Isend(B_remote, nl, MPI_BODY, left, 99, comm1d, &request[0]);
            if (save && mype == 0)
            {
                for (int i = 0; i < nl; ++i)
                    fprintf(step_file, "%f %f %f\n", B_remote[i].x, B_remote[i].y, B_remote[i].z);
            }
            MPI_Irecv(B_remote, nl, MPI_BODY, right, 99, comm1d, &request[1]);
            MPI_Waitall(2, request, status);
        }

        if (save && mype == 0)
            fclose(step_file);

        MPI_Barrier(comm1d);
        for (int i = 0; i < nl; i++)
        { /* compute new position */
            B_local[i].x += B_local[i].vx * dt;
            B_local[i].y += B_local[i].vy * dt;
            B_local[i].z += B_local[i].vz * dt;
            // bounce back without speed loss
            if (B_local[i].x > L || B_local[i].x < -L)
            {
                B_local[i].vx *= -1;
            }

            if (B_local[i].y > L || B_local[i].y < -L)
            {
                B_local[i].vy *= -1;
            }
            if (B_local[i].z > L || B_local[i].z < -L)
            {
                B_local[i].vz *= -1;
            }
        }
    }
    MPI_Barrier(comm1d);
    if (mype == 0)
    {
        double totalTime = GetTimer() / 1000.0;
        double avgTime = totalTime / (double)(nIters);
        printf("M = [%d]; avgTime: %f; totTime: %f\n", M, avgTime, totalTime);
    }
    free(B_local);
    free(B_remote);
    MPI_Finalize();
}