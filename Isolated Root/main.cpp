#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <conio.h>
#include <limits>

using namespace std;

// Исследуемая функция
double f(double x) {
    return sin(exp(x / 2.0)) - exp(-x) + 1;
}

// Производная исследуемой функции
double df(double x) {
    return (1.0 / 2.0) * cos(exp(x / 2.0)) * exp(x / 2.0) + exp(-x);
}

// Метод дихотомии
void dichotomy(double a, double b, double eps, int &n) {
    double x1;
    n++;
    if (fabs(b - a) > eps) {
        x1 = (a + b) / 2.0;
        if (f(a) * f(x1) <= 0) {
            dichotomy(a, x1, eps, n);
        } else if (f(b) * f(x1) < 0) {
            dichotomy(x1, b, eps, n);
        }
    } else if (fabs(b - a) <= eps)
        cout << "Dichotomy method root of the equation x = " << (a + b) / 2.0 << ", for E = " << eps << ", iterations n = " << n
             << endl;

}

// Метод секущих
void secant(double a, double b, double eps, int &n) {
    double k, c;
    n++;
    if (fabs(b - a) > eps) {
        k = (f(b) - f(a)) / (b - a);
        c = a - f(a) / k;
        if (f(a) * f(c) <= 0 && (f(b) * f(c) <= 0)) {
            if (fabs(a - c) < fabs(b - c)) secant(a, c, eps, n);
            else secant(c, b, eps, n);
        } else if (f(a) * f(c) <= 0) {
            secant(a, c, eps, n);
        } else if (f(b) * f(c) < 0) {
            secant(c, b, eps, n);
        }
    } else if (fabs(b - a) <= eps)
        cout << "Secant method root of the equation x = " << (a + b) / 2.0 << ", for E = " << eps << ", iterations n = " << n
             << endl;
}

// Метод квадратичной интрополяции
void quadraticInterpl(double x_, double _x, double eps, int n) {
    if (fabs(x_ - _x) > eps) {
        n++;
        double h = fabs(x_ - _x) / 2.0;
        double x0 = (x_ + _x) / 2.0;
        double a = (f(x_) - 2.0 * f(x0) + f(_x)) / (2.0 * pow(h, 2.0));
        double b = (f(x_) - f(_x)) / (2.0 * h);
        double c = f(x0);
        double d = pow(b, 2.0) - 4.0 * a * c;
        double x1 = x0 + (-b + sqrt(d)) / (2.0 * a);
        double x2 = x0 + (-b - sqrt(d)) / (2.0 * a);
        double x;
        if (x_ < x2 && x2 < _x) {
            x = x2;
            if (f(x_) * f(x) <= 0) {
                quadraticInterpl(x_, x, eps, n);
            } else if (f(x) * f(_x) < 0) {
                quadraticInterpl(x, _x, eps, n);
            }
        }
        if (x_ < x1 && x1 < _x) {
            x = x1;
            if (f(x_) * f(x) <= 0) {
                quadraticInterpl(x_, x, eps, n);
            } else if (f(x) * f(_x) < 0) {
                quadraticInterpl(x, _x, eps, n);
            }
        }
    } else if (fabs(x_ - _x) <= eps)
        cout << "Quadratic interpolation method root of the equation x = " << (x_ + _x) / 2.0 << ", for E = " << eps
             << ", iterations n = " << n
             << endl;
}

// Метод Ньютона
void newton(double x, double eps, int n)
{
    double func = f(x), dFunc;
    while (fabs(func) > eps)
    {
        class Zero{};
        func = f(x);
        dFunc = df(x);
        try {
            if (dFunc==0) throw Zero();
            x = x - func / dFunc;
            n++;
        }
        catch (const Zero & ex) {
            cout << "Division by zero\n";
            exit(1);
        }
    }
    cout << "Newton method root of the equation x = " << x << ", for E = " << eps << ", iterations n = " << n << endl;
}

// Метод последовательных приближений
void successiveApprox(double x, double eps, int n)
{
    double func = numeric_limits<double>::max();
    while (fabs(func - x) > eps)
    {
        func = x;
        x = f(x);
        n++;
    }
    cout << "Successive approximations method root of the x = f(x) equation x = " << x << ", for E = " << eps << ", iterations n = " << n << endl;
}

// Тело функции
int main() {
    double xLeft, xRight, epsilon = 0.00000001 /*numeric_limits<double>::epsilon()*/;
    int n;
    cout << "Enter function endpoints" << endl;
    cin >> xLeft >> xRight;
    // 2.5 and 3
    // Настройки точности вывода
    const auto digits = 8;
    cout << setw(digits);
    cout << fixed << std::setprecision(digits);
    // Применение методов
    n = 0;
    dichotomy(xLeft, xRight, epsilon, n);
    n = 0;
    secant(xLeft, xRight, epsilon, n);
    n = 0;
    quadraticInterpl(xLeft, xRight, epsilon, n);
    n = 0;
    newton((xLeft + xRight) / 2.0, epsilon, n);
    n = 0;
    successiveApprox(xLeft, epsilon, n);
}