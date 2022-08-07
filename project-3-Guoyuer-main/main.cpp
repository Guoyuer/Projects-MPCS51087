#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
using namespace std;
struct Vector {
  // cannot use STL in CUDA, so...
  double data[3];
  Vector() {
    for (int i = 0; i < 3; i++)
      data[i] = 0;
  };
  Vector(double x, double y, double z) {
    data[0] = x;
    data[1] = y;
    data[2] = z;
  }
  Vector(const Vector &v) {
    for (int i = 0; i < 3; i++)
      data[i] = v.data[i];
  };
  double &operator[](int idx) { return data[idx]; }
  Vector scale(double s) {
    Vector res;
    for (int i = 0; i < 3; i++)
      res[i] = data[i] * s;
    return res;
  }

  double length() {
    double sum = 0;
    for (int i = 0; i < 3; i++)
      sum += data[i] * data[i];
    return sqrt(sum);
  }

  Vector normalize() {
    Vector res = *this;
    double len = this->length();
    return res.scale(1 / len);
  }

  double operator*(Vector &v) {
    double sum = 0;
    for (int i = 0; i < 3; i++) {
      sum += data[i] * v[i];
    }
    return sum;
  }
  friend ostream &operator<<(ostream &os, const Vector &v);
};

ostream &operator<<(ostream &os, const Vector &v) {
  os << "(" << v.data[0] << ", " << v.data[1] << ", " << v.data[2] << ')'
     << endl;
  return os;
}
Vector operator-(Vector &a, Vector &b) {
  Vector res;
  for (int i = 0; i < 3; i++) {
    res[i] = a[i] - b[i];
  }
  return res;
}


bool check(Vector &W, Vector &V, Vector &C, double Wmax, double R,
           double second) {

  return abs(W[0]) < Wmax && abs(W[2]) < Wmax && second > 0;
}
void save_matrix(double **A, int N, string &fname) {
  FILE *fp = fopen(fname.c_str(), "w");

  for (long i = 0; i < N; i++) {
    for (long j = 0; j < N; j++) {
      fprintf(fp, "%.9le ", A[i][j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
}

int main(int argc, char *argv[]) {
  Vector L = Vector(4, 4, -1);
  Vector C = Vector(0, 12, 0);
  double R = 6, Wy = 10, Wmax = 10;
  int n_ray = stoi(argv[1]);
  int n_grid_1d = stoi(argv[2]);

  std::random_device rd;
  std::mt19937 e2(rd());
  std::uniform_real_distribution<> dist1(0, 2 * M_PI);
  std::uniform_real_distribution<> dist2(-1, 1);

  double **G = (double **)malloc(n_grid_1d * sizeof(double *));
  for (int i = 0; i < n_grid_1d; i++) {
    G[i] = (double *)malloc(n_grid_1d * sizeof(double));
    memset(G[i], 0, sizeof(double) * n_grid_1d);
  }
  Vector I, N, S;
  for (int n = 1; n <= n_ray; n++) {
    Vector W, V;
    double second;
    do {
      double phi = dist1(e2);
      double cos_th = dist2(e2);
      double sin_th = sqrt(1 - cos_th * cos_th);
      double Vx = sin_th * cos(phi);
      double Vy = sin_th * sin(phi);
      double Vz = cos_th;
      V = Vector(Vx, Vy, Vz);
      W = V.scale(Wy / Vy);
      second = pow((V * C), 2) + pow(R, 2) - C * C;
    } while (!check(W, V, C, Wmax, R, second));
    // cout << W <<endl;
    double t = (V * C) - sqrt(second);
    I = V.scale(t);
    N = (I - C).normalize();
    S = (L - I).normalize();
    double b = max(0.0, S * N);
    int i = W[0] / (Wmax / (n_grid_1d / 2)) + n_grid_1d / 2;
    int j = W[2] / (Wmax / (n_grid_1d / 2)) + n_grid_1d / 2;
    G[n_grid_1d - 1 - i][j] += b;
  }
  string filename("window.out");
  save_matrix(G, n_grid_1d, filename);

  return 0;
}