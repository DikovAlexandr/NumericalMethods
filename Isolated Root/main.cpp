#include <iostream>
#include <cmath>
#include <conio.h>
#include <limits>

// рассмотреть случай когда корней нет и когда корень поадает на границу

using namespace std;

double f(double x) {
    return sin(exp(x / 2)) - exp(-x);
}

void dichotomy(double a, double b, double eps, int &n) {
    double x1;
    n++;
    //cout << a << " " << b << " " << (fabs(b - a) > eps) << endl;
    if (fabs(b - a) > eps) {
        x1 = (a + b) / 2;
        if (f(a) * f(x1) <= 0) {
            dichotomy(a, x1, eps, n);
        } else if (f(b) * f(x1) < 0) {
            dichotomy(x1, b, eps, n);
        }
    } else if (fabs(b - a) <= eps)
        cout << "Dichotomy root of the equation x = " << (a + b) / 2 << ", for E = " << eps << ", iterations n = " << n
             << endl;
}

void secant(double a, double b, double eps, int &n) {
    double k, c;
    n++;
    //cout << a << " " << b << " " << (fabs(b - a) > eps) << endl;
    if (fabs(b - a) > eps) {
        k = (f(b) - f(a)) / (b - a);
        c = a - f(a) / k;
        if (f(a) * f(c) <= 0 && (f(b) * f(c) <= 0)) {
            if (fabs(a - c) < fabs(b - c)) secant(a, c, eps, n);
            else secant(c, b, eps, n);
        } else if (f(a) * f(c) < 0) {
            secant(a, c, eps, n);
        } else if (f(b) * f(c) < 0) {
            secant(c, b, eps, n);
        }
    } else if (fabs(b - a) <= eps)
        cout << "Secant root of the equation x = " << (a + b) / 2 << ", for E = " << eps << ", iterations n = " << n
             << endl;
}

void quadraticInterpl(double x_, double _x, double eps, int n) {
    if (fabs(x_ - _x) > eps) {
        n++;
        double h = fabs(x_ - _x) / 2;
        double x0 = (x_ + _x) / 2;
        double a = (f(x_) - 2 * f(x0) + f(_x)) / (2 * pow(h, 2));
        double b = (f(x_) - f(_x)) / (2 * h);
        double c = f(x0);
        double d = pow(b, 2) - 4 * a * c;
        double x1 = x0 + (-b + sqrt(d)) / (2 * a);
        double x2 = x0 + (-b - sqrt(d)) / (2 * a);
        double x = 0;
        //cout << x_ << " " << _x << " " << endl;
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
        cout << "Quadratic Interpolation root of the equation x = " << (x_ + _x) / 2 << ", for E = " << eps
             << ", iterations n = " << n
             << endl;
}

int main() {
    double xLeft, xRight, epsilon = 0.0000001 /*numeric_limits<double>::epsilon()*/, h;
    int n;
    cout << "Enter endpoints" << endl;
    cin >> xLeft >> xRight;
    n = 0;
    dichotomy(xLeft, xRight, epsilon, n);
    n = 0;
    secant(xLeft, xRight, epsilon, n);
    n = 0;
    quadraticInterpl(xLeft, xRight, epsilon, n);
}