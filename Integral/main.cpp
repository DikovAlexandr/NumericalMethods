#include <iostream>
#include <cmath>

using namespace std;

double f(double x) {
    //return (pow(x, 4) + 5 * pow(x, 3) + 2 * pow(x, 2) + 15 * x + 2);
    return exp(-x);
}

double analyticIntegral(double a, double b) {
    //return pow(b, 5) / 5 + 5 * pow(b, 4) / 4 + 2 * pow(b, 3) / 3 + 15 * pow(b, 2) / 2 + 2 * b;
    return 1 - exp(-3);
}

double midpointSum(double a, double b, double h) {
    double integral;
    double n; // задаём число разбиений n
    n = (b - a) / h;
    // вычисляем интеграл по формуле центральных прямугольников
    integral = 0.0;
    for (int i = 1; i <= n; i++)
        integral = integral + h * f(a + h * (i - 0.5));
    return integral;
}

double trapezoidalSum(double a, double b, double h) {
    double integral;
    double n; // задаём число разбиений n
    n = (b - a) / h;
    // вычисляем интеграл по формуле трапеций
    integral = h * (f(a) + f(b)) / 2.0;
    for (int i = 1; i <= n - 1; i++)
        integral = integral + h * f(a + h * i);
    return integral;
}

double simpsonSum(double a, double b, double h) {
    double integral;
    double n; // задаём число разбиений n
    n = (b - a) / h;
    // вычисляем интеграл по формуле Симпсона
    integral = h * (f(a) + f(b)) / 6.0;
    for (int i = 1; i <= n; i++)
        integral = integral + 4.0 / 6.0 * h * f(a + h * (i - 0.5));
    for (int i = 1; i <= n - 1; i++)
        integral = integral + 2.0 / 6.0 * h * f(a + h * i);
    return integral;
}

double threeEighthSum(double a, double b, double h) {
    double n = (b - a) / h;
    double integral = (f(0) - f(h)) * 3 * h / 8;
    for (int i = 1; i <= n; i++) {
        if (i % 3 != 0)
            integral += 3 * 3 * h / 8 * f(h * i);
        else integral += 2 * 3 * h / 8 * f(h * i);
    }
    return integral;
}

double midpointAccuracy(double h) {
    return pow(h, 3) / 6.0;
}

/*double trapezoidalAccuracy(double h) {
    return h*h/12*f''
}*/

double threeEighthsAccuracy(double h) {
    return pow(h, 5) * 24 * 3 / 80;
}

double simpsonAccuracy(double h) {
    return pow(h, 5) * 24 / 90;
    // где 24 это четвертая производная от f
}

int main() {
    double a = 0;
    double b = 3.0;
    double h = 0.001;
    cout << "Midpoint sum " << midpointSum(a, b, h) << " different from analytic  "
         << fabs(analyticIntegral(a, b) - midpointSum(a, b, h)) << " while the remainder is "
         << midpointAccuracy(h) << endl;
    cout << "Trapezoidal sum " << trapezoidalSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - trapezoidalSum(a, b, h)) << endl;
    cout << "Three eighths sum " << threeEighthSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - threeEighthSum(a, b, h)) << " while the remainder is "
         << threeEighthsAccuracy(h) << " bool "
         << (fabs(analyticIntegral(a, b) - threeEighthSum(a, b, h)) < threeEighthsAccuracy(h)) << endl;
    cout << "Simpson sum " << simpsonSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - simpsonSum(a, b, h)) << " while the remainder is "
         << simpsonAccuracy(h) << " bool " << (fabs(analyticIntegral(a, b) - simpsonSum(a, b, h)) < simpsonAccuracy(h))
         << endl;
    cout << "Analytic result " << analyticIntegral(a, b) << endl;
    return 0;
}
