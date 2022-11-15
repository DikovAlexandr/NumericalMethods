#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Решаем уравнение y' = f(x,y)
double f(double x, double y) {
    return -1.0 * y;
}

// Аналитическое решение y = epx(-x)
double y(double x) {
    return exp(-x);
}

// Метод Рунге-Кутты 2-го порядка
void RK2(double xn, double x, double y0, double h, double w) {
    ofstream fout;
    fout.open(R"(..\\..\\NumericalMethods\\RungeKutta\\RK2.csv)");
    int n = ceil((x - xn) / h);
    double k1, k2, k3, k4, yn = y0;
    for (int i = 1; i <= n; i++) {
        fout << xn << " " << yn << " " << y(xn) << endl;
        k1 = f(xn, yn);
        k2 = f(xn + h / (2 * w) , yn + k1 / (2 * w));
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

// Поскольку численное решение не позволяет найти аналитический вид функции y(x) результатом является таблица значений
// функции на заданном интервале, также должна быть задана начальная точка y(x0) = y0
int main() {
    // Задаем границы интервала значений по x, начальное условие, шаг и w
    double x0 = 0, x = 2, y0 = 1, h = 0.01, w = 0.5;
    //cout << "Runge-Kutta 4th Order method" << endl;
    RK2(x0, x, y0, h, w);
    RK4(x0, x, y0, h);
    return 0;
}