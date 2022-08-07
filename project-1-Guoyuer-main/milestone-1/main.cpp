#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>

using namespace std;

void save_matrix(vector<vector<double>> &matrix, int N, const string &filename);

int main(int argc, char **argv) {
    if (argc < 7) {
        cout << "need 6 arguments" << endl;
        return 0;
    }
    int N, NT;
    double L, T, u, v;
    N = stoi(argv[1]);
    cout << "N = " << N << endl;

    NT = stoi(argv[2]);
    cout << "NT = " << NT << endl;

    L = atof(argv[3]);
    cout << "L = " << L << endl;

    T = atof(argv[4]);
    cout << "T = " << T << endl;

    u = atof(argv[5]);
    cout << "u = " << u << endl;


    v = atof(argv[6]);
    cout << "v = " << v << endl;

    /* a double uses 8 bytes*/
    cout << "Estimated Memory Usage: " << N * N * 2 * 8 / 1024 / 1024 << "MB" << endl;
    // (1~N, 1~N) is actual cell;
    double dx = L / N, dy = L / N;
    double dt = T / NT;
    if (dt > dx / sqrt(2 * (u * u + v * v))) {
        cout << "Courant stability condition not satisfied!\n";
        return 0;
    }
    vector<vector<double>> Cn(N + 2, vector<double>(N + 2));//C(n)
    vector<vector<double>> Cnn(N + 2, vector<double>(N + 2));//C(n+1)
    /* init start*/
    double denom = 2 * (L / 4) * (L / 4);
    double x0 = L / 2, y0 = L / 2;
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            double x = (i - 1) * dx;
            double y = (j - 1) * dy;
            Cn[i][j] = exp(-(pow(x - x0, 2) + pow(y - y0, 2)) / denom);
        }
    }
    /* init end*/
    save_matrix(Cn, N, "init_matrix.csv");
    clock_t start = clock();
    for (int n = 1; n <= NT; n++) {

        /* create periodic boundaries */
        for (int k = 1; k <= N; k++) {
            Cn[k][0] = Cn[k][N];//left
            Cn[k][N + 1] = Cn[k][1];//right
            Cn[N + 1][k] = Cn[1][k];//bottom
            Cn[0][k] = Cn[N][k];//top
        }
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                Cnn[i][j] = (Cn[i - 1][j] + Cn[i + 1][j] + Cn[i][j - 1] + Cn[i][j + 1]) / 4 -
                            dt * (u * (Cn[i + 1][j] - Cn[i - 1][j]) + v * (Cn[i][j + 1] - Cn[i][j - 1])) / (2 * dx);
            }
        }
//        if (n == NT / 2) {
//            save_matrix(Cn, N, "halfway_matrix.csv");
//        }
//        if (n == NT) {
//            save_matrix(Cn, N, "final_matrix.csv");
//        }
        swap(Cnn, Cn);
    }
    clock_t end = clock();
    double time_elapsed = (double) (end - start) / CLOCKS_PER_SEC;
    printf("Time taken: %.2fs\n", time_elapsed);
    printf("Grind Rate: %d step/s\n", (int) (NT / time_elapsed));

    return 0;
}

void save_matrix(vector<vector<double>> &matrix, int N, const string &filename) {
    ofstream file;
    file.open(filename);
    file << scientific;
    for (int i = 1; i <= N; i++) {
        for (int j = 1; j <= N; j++) {
            file << matrix[i][j] << ", ";
        }
        file << "\n";
    }
    file.close();
}