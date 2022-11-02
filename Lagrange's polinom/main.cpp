#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Исследуемая функция
double f(double x) {
    return sin(1 / (x + 1));
}

// Получение значений L(x)
double lagrange(const double *x, const double *y, int n, double _x) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        double l = 1.0;
        for (int j = 0; j < n; j++)
            if (j != i) {
                l *= (_x - x[j]) / (x[i] - x[j]);
            }
        result += l * y[i];
    }
    return result;
}

// Тело функции
int main() {
    ofstream fout;
    fout.open(R"(..\\..\\NumericalMethods\\Lagrange's polinom\\plot.csv)");
    int n = 10;
    int N = n * 100;
    double x[n], y[n], L[N];
    for (int i = 0; i < n; i++) {
        x[i] = i;
        y[i] = f(x[i]);
    }
    for (int i = 0; i < N; i++) {
        double t = i / 100.0;
        L[i] = lagrange(x, y, n, t);
        fout << t << " " << L[i] << " " << f(t) << endl;
    }
    fout.close();
    return 0;
}