#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Собственная функция в бесконечноглубокой прямоугольной потенцальной яме
double psiK(double x, double a, int k) {
    return sqrt(2 / a) * sin(M_PI * k * x / a);
}

// Собственная функция в бесконечноглубокой прямоугольной потенцальной яме
double psiM(double x, double a, int m) {
    return sqrt(2 / a) * sin(M_PI * m * x / a);
}

// Энергия n-го состония частицы в бесконечноглубокой прямоугольной потенцальной яме
double Energy(double plank, double a, int n, double M) {
    return pow(M_PI, 2) * pow(plank, 2) * pow(n, 2) / (2 * M * pow(a, 2));
}

// Элемент матрицы с индексами k и m
double Vmk(int k, int m, int alpha, double a, double h) {
    double z = 0;
    int N = static_cast<int>((a + h) / h);
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

// Туло функции
int main() {
    // Объявим константы
    double plank = 1;
    double M = 1;
    double a = 20;
    int alpha = 1;
    double h = 0.005;
    int N = 10;

    // Создаем матрицу и заполняем ее
    Eigen::MatrixXd matrix(N, N);
    for (int k = 0; k < N; k++) {
        for (int m = 0; m < N; m++) {
            if (k == m) matrix(k, m) = Energy(plank, a, k + 1, M) + Vmk(k + 1, m + 1, alpha, a, h);
            else matrix(k, m) = Vmk(k + 1, m + 1, alpha, a, h);
        }
    }

    // Выводим матрицу для проверки
    //cout << matrix << endl;

    // Создаем класс из матрицы для нахождения ее собственный значений и собственных векторов
    EigenSolver<MatrixXd> es(matrix);

    // Поскольку собственные значения могут быть комплекснозначны будем хранить их в двух представлениях
    // Комплексном
    complex<double> complexValues[N];
    // И действительном
    double doubleValues[N];
    for (int i = 0; i < es.eigenvalues().size(); ++i) {
        complexValues[i] = es.eigenvalues()[i];
        doubleValues[i] = complexValues[i].real();
        // cout << doubleValues[i] << endl;
    }

    // cout << es.eigenvectors().cols() << endl;

    double doubleVectors[N][N];
    for (int i = 0; i < es.eigenvectors().cols(); ++i) {
        for (int j = 0; j < N; ++j) {
            doubleVectors[i][j] = es.eigenvectors().col(i)[j].real() * doubleValues[i];
            cout << doubleVectors[i][j] << " ";
        }
        cout << endl << endl;
    }

    // Получим значения волновой функции
    int aN = static_cast<int>((a + h )/ h);
    double x[aN];
    for (int i = 0; i < aN; ++i) {
        x[i] = h * i;
    }

    ofstream fout;
    fout.open(R"(..\\..\\Matrix Schrodinger\\plot.csv)");

    // Массив значений волновой функции
    double psi[aN];
    // Нужный номер состояния
    int n = 1;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < aN; ++j) {
            psi[j] += psiK(x[j], a, i) * doubleVectors[n][i];
            // cout << psi[j] << endl;
        }
    }
    for (int j = 0; j < aN; ++j) {
        fout << x[j] << " " << psi[j] << endl;
    }
    fout.close();
}
