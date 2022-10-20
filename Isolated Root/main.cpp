#include <iostream>
#include <cmath>
#include <conio.h>
#include <limits>

// рассмотреть случай когда корней нет и когда корень поадает на границу

using namespace std;

double f(double x) {
    return sin(x / 2);
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
        c = (k * a - f(a)) / k;
        if (f(a) * f(c) < 0 && (f(b) * f(c) < 0)) {
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

void quadraticInterpl(double x_, double _x) {
    double h = fabs(x_ - _x) / 2;
    double x0 = (x_ + _x) / 2;
    double a = (f(x_) - 2 * f(x0) + f(_x)) / (2 * pow(h, 2));
    double b = (f(x_) - f(_x)) / (2 * h);
    double c = f(x0);
    double d = pow(b, 2) - 4 * a * c;
    cout << a << " " << b << " " << c << endl;
    cout << "Quadratic interpolation root of the equation x = " << (-b + sqrt(d)) / (2 * a) << ", for h = " << h
         << endl;
}

int main() {
    double xLeft, xRight, epsilon = 0.0001 /*numeric_limits<double>::epsilon()*/, h;
    int n;
    cout << "Enter endpoints" << endl;
    cin >> xLeft >> xRight;
    n = 0;
    dichotomy(xLeft, xRight, epsilon, n);
    n = 0;
    secant(xLeft, xRight, epsilon, n);
    h = 0.01;
    n = 0;
    quadraticInterpl(xLeft, xRight);
}