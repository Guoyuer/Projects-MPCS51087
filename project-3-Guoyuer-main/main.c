#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#define M_PI 3.14159265358979323846 /* pi */
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define NTHREADS_PER_BLOCK 256
typedef struct vector
{
  double x;
  double y;
  double z;
} vector;

vector scale(vector v, double s)
{
  vector nv;
  nv.x = v.x * s;
  nv.y = v.y * s;
  nv.z = v.z * s;
  return nv;
}

vector subtract(vector a, vector b)
{
  vector nv;
  nv.x = a.x - b.x;
  nv.y = a.y - b.y;
  nv.z = a.z - b.z;
  return nv;
}

double product(vector a, vector b)
{
  double sum = 0;
  sum += a.x * b.x;
  sum += a.y * b.y;
  sum += a.z * b.z;
  return sum;
}

double length(vector v)
{
  return sqrt(product(v, v));
}
vector normalize(vector v)
{
  double len = length(v);
  return scale(v, 1.0 / len);
}

void print(vector *v)
{
  printf("(%.2f, %.2f, %.2f)\n", v->x, v->y, v->z);
}

int check(vector W, double Wmax, double R, double term)
{
  return (abs(W.x) < Wmax && abs(W.z) < Wmax && term > 0);
}
void save_matrix(double **A, int N, char *fname)
{
  FILE *fp = fopen(fname, "w");

  for (long i = 0; i < N; i++)
  {
    for (long j = 0; j < N; j++)
    {
      fprintf(fp, "%.9le ", A[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}
double randfrom(double min, double max)
{
  double range = (max - min);
  double div = RAND_MAX / range;
  return min + (rand() / div);
}
int main(int argc, char *argv[])
{
  // vector L = create_vec(4, 4, -1);
  vector L = {4, 4, -1};
  vector C = {0, 12, 0};
  double R = 6, Wy = 10, Wmax = 10;
  int n_ray = atoi(argv[1]);
  int n_grid_1d = atoi(argv[2]);

  double **G = (double **)malloc(n_grid_1d * sizeof(double *));
  for (int i = 0; i < n_grid_1d; i++)
  {
    G[i] = (double *)malloc(n_grid_1d * sizeof(double));
    memset(G[i], 0, sizeof(double) * n_grid_1d);
  }

  vector I, N, S;
  time_t t;
  srand((unsigned)time(&t));

  clock_t start, end;
  start = clock();

  for (int n = 1; n <= n_ray; n++)
  {
    vector W, V;
    double term;
    do
    {
      double phi = randfrom(0, 2 * M_PI);
      double cos_th = randfrom(-1, 1);
      double sin_th = sqrt(1 - cos_th * cos_th);
      double Vx = sin_th * cos(phi);
      double Vy = sin_th * sin(phi);
      double Vz = cos_th;
      V.x = Vx;
      V.y = Vy;
      V.z = Vz;
      W = scale(V, Wy / Vy);
      term = pow(product(V, C), 2) + pow(R, 2) - product(C, C);
    } while (!check(W, Wmax, R, term));

    double t = product(V, C) - sqrt(term);
    I = scale(V, t);
    vector I_C = subtract(I, C);
    N = normalize(I_C);
    vector L_I = subtract(L, I);
    S = normalize(L_I);
    double b = MAX(0.0, product(S, N));
    int i = W.x / (Wmax / (n_grid_1d / 2)) + n_grid_1d / 2;
    int j = W.z / (Wmax / (n_grid_1d / 2)) + n_grid_1d / 2;
    G[n_grid_1d - 1 - i][j] += b;
  }
  end = clock();
  double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC; // in seconds
  printf("time elapsed on CPU: %.6f(s)", time_taken);
  save_matrix(G, n_grid_1d, "window.out");
  return 0;
}