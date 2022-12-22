#include <iostream>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Собственная функция в бесконечноглубокой прямоугольной потенцальной яме
double psiK(double x, double a, int k) {
    return sqrt(2 / a) * sin(M_PI * k * x / a);
}

double psiM(double x, double a, int m) {
    return sqrt(2 / a) * sin(M_PI * m * x / a);
}

double Energy(double plank, double a, int n, double M) {
    return pow(M_PI, 2) * pow(plank, 2) * pow(n, 2) / (2 * M * pow(a, 2));
}

double Vmk(int k, int m, int alpha, double a, double h) {
    double z = 0;
    int N = static_cast<int>(a + h / h);
    double x[N];
    for (int i = 0; i < N; ++i) {
        x[i] = h * i;
    }
    double y[N];
    for (int i = 0; i < N; ++i) {
        y[i] = pow(x[i], 2) * alpha * psiK(x[i], a, k) * psiM(x[i], a, m);
    }
    for (int i = 0; i < N - 1; ++i) {
        z += h * (y[i] + y[i + 1]) / 2;
    }
    return z;
}

int main() {
    double plank = 1;
    double M= 1;
    double a = 10;
    int alpha = 1;
    double h = 0.005;
    int N=10;
    Eigen::MatrixXd matrix(N,N);
    for (int k = 0; k < N; ++k) {
        for (int m = 0; m < N; ++m) {
            if (k==m) matrix(k, m) = Energy(plank, a, k, M) + Vmk(k, m, alpha, a, h);
            else matrix(k, m) = Vmk(k, m, alpha, a, h);
        }
    }
    matrix.eigenvalues();
}
