#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "timer.h"
#include <time.h> //for random seed
#include <omp.h>

#define SOFTENING 1e-3f
#define G 9.8f
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
void init(Particle *p, int n)
{
    for (int i = 0; i < n; i++)
    {
        p[i].m = 1.0f;

        p[i].x = getRandIn(-5, 5);
        p[i].y = getRandIn(-2, 2);
        p[i].z = 0.0f; // on same plane: Z = 0
        p[i].vx = 0.0f;
        p[i].vy = 0.0f;
        p[i].vz = 0.0f;
    }
}

/* calculate all interparticle forces and update instantaneous velocities */
void calc_force(Particle *p, float dt, int n)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        double Fx = 0.0f;
        double Fy = 0.0f;
        double Fz = 0.0f;

        for (int j = 0; j < n; j++)
        {
            if (i == j)
                continue;
            /* calculate net particle for on i'th particle */
            float dx = p[j].x - p[i].x;
            float dy = p[j].y - p[i].y;
            float dz = p[j].z - p[i].z;
            float D_rij = sqrtf(dx * dx + dy * dy + dz * dz + SOFTENING); // distance of two bodies
            float D_rij3 = D_rij * D_rij * D_rij;
            float left_term = G * p[i].m * p[j].m / D_rij3;

            Fx += dx * left_term;
            Fy += dy * left_term;
            Fz += dz * left_term;
        }
        /* update instantaneous velocity based on force and timestep */
        p[i].vx += dt * Fx;
        p[i].vy += dt * Fy;
        p[i].vz += dt * Fz;
    }
}

int main(const int argc, const char **argv)
{
    srand(42);
    int nParticles = atoi(argv[1]); /* number of particles */
    int nIters = atoi(argv[2]);     /* number of steps in simulation */
    int save = atoi(argv[3]);       /* if save == 1 then output to file*/
    double L = 10;
    const float dt = 0.005f; /* time step   */

    Particle *p = malloc(nParticles * sizeof(Particle));

    init(p, nParticles); /* Init mass, pos and vel data */
    FILE *datafile = NULL;
    if (save)
    {
        datafile = fopen("particles_ref.dat", "w");
        fprintf(datafile, "%d %d %d\n", nParticles, nIters, 0);
    }

    /* ------------------------------*/
    /*     MAIN LOOP                 */
    /* ------------------------------*/
    StartTimer();
    for (int iter = 1; iter <= nIters; iter++)
    {
        if (save)
        {
            printf("iteration:%d\n", iter);
            for (int i = 0; i < nParticles; ++i)
                fprintf(datafile, "%f %f %f\n", p[i].x, p[i].y, p[i].z);
            fprintf(datafile, "\n");
        }

        calc_force(p, dt, nParticles); /* compute interparticle forces and update vel */

        for (int i = 0; i < nParticles; i++)
        { /* compute new position */
            p[i].x += p[i].vx * dt;
            p[i].y += p[i].vy * dt;
            p[i].z += p[i].vz * dt;
            // bounce back without speed loss
            if (p[i].x > L || p[i].x < -L)
            {
                p[i].vx *= -1;
            }

            if (p[i].y > L || p[i].y < -L)
            {
                p[i].vy *= -1;
            }
            if (p[i].z > L || p[i].z < -L)
            {
                p[i].vz *= -1;
            }
        }
    }

    double totalTime = GetTimer() / 1000.0;
    double avgTime = totalTime / (double)(nIters);

    if (save)
    {
        fclose(datafile);
    }

    printf("avgTime: %f   totTime: %f \n", avgTime, totalTime);
    free(p);
}