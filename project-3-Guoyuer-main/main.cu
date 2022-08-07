#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <cuda.h>
#include <stdint.h>
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#define M_PI 3.14159265358979323846 /* pi */
#define MAX_BLOCKS_PER_DIM 2147483647
#define MAX_RAND_NUM 200


typedef struct vector
{
  double x;
  double y;
  double z;
} vector;

__device__ uint64_t fast_forward_LCG(uint64_t seed, uint64_t n);
__device__ double LCG_random_double(uint64_t *seed);

__device__ vector scale(vector v, double s)
{
  vector nv;
  nv.x = v.x * s;
  nv.y = v.y * s;
  nv.z = v.z * s;
  return nv;
}

__device__ vector subtract(vector a, vector b)
{
  vector nv;
  nv.x = a.x - b.x;
  nv.y = a.y - b.y;
  nv.z = a.z - b.z;
  return nv;
}

__device__ double product(vector a, vector b)
{
  double sum = 0;
  sum += a.x * b.x;
  sum += a.y * b.y;
  sum += a.z * b.z;
  return sum;
}

__device__ double length(vector v)
{
  return sqrt(product(v, v));
}
__device__ vector normalize(vector v)
{
  double len = length(v);
  return scale(v, 1.0 / len);
}

void print(vector v)
{
  printf("(%.2f, %.2f, %.2f)\n", v.x, v.y, v.z);
}

__global__ void fill_vec(vector *v, double x, double y, double z)
{
  v->x = x;
  v->y = y;
  v->z = z;
}
__device__ int check(vector W, double Wmax, double R, double term);

__global__ void compute(double *G, int n_ray, int n_grid_1d, vector L, vector C,
                        double R, double Wy, double Wmax)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  // two seeds, for two distribution
  
  uint64_t init_seed1 = (uint64_t)tid * 3153 + 451028;
  uint64_t seed1 = fast_forward_LCG(init_seed1, tid * MAX_RAND_NUM);

  uint64_t init_seed2 = (uint64_t)tid * 784031 + 145648;
  uint64_t seed2 = fast_forward_LCG(init_seed2, tid * MAX_RAND_NUM);
  vector I, N, S;
  // grid-stride loop
  for(int n = 0; n < n_ray; n += blockDim.x * gridDim.x)
  {
    vector W, V;
    double term;
    do
    {
      double phi = LCG_random_double(&seed1) * 2 * M_PI;
      double cos_th = -1 + LCG_random_double(&seed2) * 2;
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
    int idx = (n_grid_1d - 1 - i) * n_grid_1d + j; // flatten the 2d matrix;
    atomicAdd(&G[idx], b);
  }
}

__device__ int check(vector W, double Wmax, double R, double term)
{
  return (abs(W.x) < Wmax && abs(W.z) < Wmax && term > 0);
}
void save_matrix(double *A, int N, char *fname)
{
  FILE *fp = fopen(fname, "w");

  for (long i = 0; i < N; i++)
  {
    for (long j = 0; j < N; j++)
    {
      long idx = i * N + j;
      fprintf(fp, "%.7le ", A[idx]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

__device__ uint64_t fast_forward_LCG(uint64_t seed, uint64_t n)
{
  const uint64_t m = 9223372036854775808ULL;
  uint64_t a = 2806196910506780709ULL;
  uint64_t c = 1ULL;
  n = n % m;
  uint64_t a_new = 1;
  uint64_t c_new = 0;

  while (n > 0)
  {
    if (n & 1)
    {
      a_new *= a;
      c_new = c_new * a + c;
    }
    c *= (a + 1);
    a *= a;
    n >>= 1;
  }
  return (a_new * seed + c_new) % m;
}

__device__ double LCG_random_double(uint64_t *seed)
{
  const uint64_t m = 9223372036854775808ULL;
  const uint64_t a = 2806196910506780709ULL;

  const uint64_t c = 1ULL;
  *seed = (a * (*seed) + c) % m;
  return (double)(*seed) / (double)m;
}

int main(int argc, char *argv[])
{

  int n_ray = atoi(argv[1]);
  int n_grid_1d = atoi(argv[2]);
  int NTHREADS_PER_BLOCK = atoi(argv[3]);
  if(NTHREADS_PER_BLOCK > 1024){
    printf("The maximum number of thread per block is 1024!\n");
    return 0;
  }
  vector L = {4, 4, -1};
  vector C = {0, 12, 0};
  double R = 6, Wy = 10, Wmax = 10;
  size_t len = n_grid_1d * n_grid_1d * sizeof(double);
  double *G_h = (double *)malloc(len);
  double *G;
  cudaMalloc((void **)&G, len);
  int nblocks = MIN(n_ray / NTHREADS_PER_BLOCK + 1, MAX_BLOCKS_PER_DIM);
  printf("nblocks: %d;  NTHREADS_PER_BLOCK: %d; avg_nrays_per_thread %.2f\n", nblocks, NTHREADS_PER_BLOCK, (double) n_ray/(nblocks*NTHREADS_PER_BLOCK));
  // if (nblocks == MAX_BLOCKS_PER_DIM)
  // {
  //   printf("Total number of thread < number of rays!\n");
  //   printf("Still works!\n");
  // }
  cudaDeviceSynchronize();

  cudaEvent_t start_device, stop_device;
  cudaEventCreate(&start_device);
  cudaEventCreate(&stop_device);
  float time_device;
  cudaEventRecord(start_device, 0);

  compute<<<nblocks, NTHREADS_PER_BLOCK>>>(G, n_ray, n_grid_1d, L, C, R, Wy, Wmax);

  cudaEventRecord(stop_device, 0);
  cudaEventSynchronize(stop_device);
  cudaEventElapsedTime(&time_device, start_device, stop_device);
  printf("time elapsed on GPU: %f(s)\n", time_device / 1000.);
  cudaEventDestroy(start_device);
  cudaEventDestroy(stop_device);

  cudaDeviceSynchronize();

  cudaMemcpy(G_h, G, len, cudaMemcpyDeviceToHost);

  save_matrix(G_h, n_grid_1d, "window.out");
  return 0;
}