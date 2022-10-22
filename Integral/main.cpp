#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <conio.h>
#include <limits>

using namespace std;

// Исследуемая функция
double f(double x) {
    // return (pow(x, 4) + 5 * pow(x, 3) + 2 * pow(x, 2) + 15 * x + 2);
    // return exp(-x);
    return x;
}

// Аналитический интеграл
double analyticIntegral(double a, double b) {
    // return pow(b, 5) / 5 + 5 * pow(b, 4) / 4 + 2 * pow(b, 3) / 3 + 15 * pow(b, 2) / 2 + 2 * b;
    // return 1 - exp(-3);
    return 0.5;
}

// Метод центральных прямоугольников
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

// Метод трапеций
double trapezoidalSum(double a, double b, double h) {
    double integral;
    double n; // задаём число разбиений n
    n = (b - a) / h;
    // вычисляем интеграл по формуле трапеций
    integral = h * (f(a) + f(b)) / 2.0;
    for (int i = 1; i < n; i++)
        integral = integral + h * f(a + h * i);
    return integral;
}

// Метод трех восьмых
double threeEighthSum(double a, double b, double h) {
    double n = (b - a) / h;
    double integral = (f(a) - f(b)) * 3.0 * h / 8.0;
    for (int i = 1; i < n; i++) {
        if (i % 3 != 0)
            integral += 3.0 * 3.0 * h / 8.0 * f(a + h * i);
        else integral += 2.0 * 3.0 * h / 8.0 * f(a + h * i);
    }
    return integral;
}

// Метод Симпсона
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

// Оценка остаточных членов
/*
// Теоретическая погрешность метода прямоугольников
double midpointAccuracy(double h) {
    return pow(h, 3) / 6.0;
}

// Теоретическая погрешность метода трапеций
double trapezoidalAccuracy(double h) {
    return pow(h, 2) * (b - a) /12.0 * maxF2;
    // где maxF4 это максимальная величина второй производной от f
}

// Теоретическая погрешность метода трех восьмых
double threeEighthsAccuracy(double h) {
    return pow(h, 5) * maxF4 * 3.0 / 80.0;
    // где maxF4 это максимальная величина четвертой производной от f
}

// Теоретическая погрешность метода Симпсона
double simpsonAccuracy(double h) {
    return pow(h, 4) * (b - a) * maxF4 / 180;
    // где maxF4 это максимальная величина четвертой производной от f
}
*/
int main() {
    double a = 0.0;
    double b = 1.0;
    double h = 0.001;
    /*
    // Настройки точности вывода
    const auto digits = 8;
    cout << setw(digits);
    cout << fixed << std::setprecision(digits);
    */
    // Результаты вычислений
    cout << "Midpoint sum " << midpointSum(a, b, h) << " different from analytic  "
         << fabs(analyticIntegral(a, b) - midpointSum(a, b, h)) << endl;
    cout << "Trapezoidal sum " << trapezoidalSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - trapezoidalSum(a, b, h)) << endl;
    cout << "Three eighths sum " << threeEighthSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - threeEighthSum(a, b, h)) << endl;
    cout << "Simpson sum " << simpsonSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - simpsonSum(a, b, h)) << endl;
    cout << "Analytic result " << analyticIntegral(a, b) << endl;
    return 0;
}
