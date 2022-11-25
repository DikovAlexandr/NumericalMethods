#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <array>

using namespace std;

// Решаем уравнение y' = f(x,y)
double f(double x, double y) {
    return y;
}

// Аналитическое решение y = epx(-x)
double y(double x) {
    return exp(x);
}

// Метод Рунге-Кутты 2-го порядка
void RK2(double xn, double x, double y0, double h, double w) {
    ofstream fout;
    fout.open(R"(..\\..\\NumericalMethods\\RungeKutta\\RK2.csv)");
    int n = ceil((x - xn) / h);
    double k1, k2, yn = y0;
    for (int i = 1; i <= n; i++) {
        fout << xn << " " << yn << " " << y(xn) << endl;
        k1 = f(xn, yn);
        k2 = f(xn + 0.5 * h / w, yn + 0.5 * h * k1 / w);
        yn += h * ((1 - w) * k1 + w * k2);
        xn += h;
    }
    fout << xn << " " << yn << " " << y(xn) << endl;
    fout.close();
}

// Метод Рунге-Кутты 4-го порядка
void RK4(double xn, double x, double y0, double h) {
    ofstream fout;
    fout.open(R"(..\\..\\NumericalMethods\\RungeKutta\\RK4.csv)");
    int n = ceil((x - xn) / h);
    double k1, k2, k3, k4, yn = y0;
    for (int i = 1; i <= n; i++) {
        fout << xn << " " << yn << " " << y(xn) << endl;
        k1 = f(xn, yn);
        k2 = f(xn + 0.5 * h, yn + 0.5 * h * k1);
        k3 = f(xn + 0.5 * h, yn + 0.5 * h * k2);
        k4 = f(xn + h, yn + h * k3);
        yn += h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
        xn += h;
    }
    fout << xn << " " << yn << " " << y(xn) << endl;
    fout.close();
}

void adams(double a, double b, double y0, double h) {
    ofstream fout;
    fout.open(R"(..\\..\\NumericalMethods\\RungeKutta\\Adams.csv)");
    double k1, k2, k3, k4;
    int n = ceil((b - a) / h);
    double xn[n];
    for (int i = 0; i < n; ++i) {
        xn[i] = a + i * h;
    }
    double yn[n];
    yn[0] = y0;
    for (int i = 0; i < 4; i++) {
        fout << xn[i] << " " << yn[i] << " " << y(xn[i]) << endl;
        k1 = f(xn[i], yn[i]);
        k2 = f(xn[i] + 0.5 * h, yn[i] + 0.5 * h * k1);
        k3 = f(xn[i] + 0.5 * h, yn[i] + 0.5 * h * k2);
        k4 = f(xn[i] + h, yn[i] + h * k3);
        yn[i + 1] = yn[i] + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
    }

    for (int i = 4; i < n; i++) {
        yn[i] = yn[i - 1] + h * (55 * f(xn[i - 1], yn[i - 1]) - 59 * f(xn[i - 1], yn[i - 1]) +
                                 37 * f(xn[i - 1], yn[i - 1]) - 9 * f(xn[i - 1], yn[i - 1])) / 24;
        fout << xn[i] << " " << yn[i] << " " << y(xn[i]) << endl;
    }
}

// Поскольку численное решение не позволяет найти аналитический вид функции y(x) результатом является таблица значений
// функции на заданном интервале, также должна быть задана начальная точка y(x0) = y0
int main() {
    // Задаем границы интервала значений по x, начальное условие, шаг и w
    double x0 = 0, x = 10, y0 = 1, h = 0.1, w = 0.5;
    RK2(x0, x, y0, h, w);
    RK4(x0, x, y0, h);
    adams(x0, x, y0, h);
    return 0;
}