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
    return sin(x);
}

// Вторая производная
double d2_f(double x) {
    // return (pow(x, 4) + 5 * pow(x, 3) + 2 * pow(x, 2) + 15 * x + 2);
    // return exp(-x);
    return -sin(x);
}

// Четвертая производная
double d4_f(double x) {
    // return (pow(x, 4) + 5 * pow(x, 3) + 2 * pow(x, 2) + 15 * x + 2);
    // return exp(-x);
    return sin(x);
}

// Аналитический интеграл
double analyticIntegral(double a, double b) {
    // return pow(b, 5) / 5 + 5 * pow(b, 4) / 4 + 2 * pow(b, 3) / 3 + 15 * pow(b, 2) / 2 + 2 * b;
    // return 1 - exp(-3);
    return 2;
}

// Максимум второй производной на промежутке интегрирования
double max_d2_f(double a, double b, double h) {
    double max = 0;
    double n; // задаём число разбиений n
    n = (b - a) / h;
    for (int i = 0; i <= n; ++i) {
        if (max < abs(d2_f(a + h * i))) max = abs(d2_f(a + h * i));
    }
    return max;
}

// Максимум четвертой производной на промежутке интегрирования
double max_d4_f(double a, double b, double h) {
    double max = 0;
    double n; // задаём число разбиений n
    n = (b - a) / h;
    for (int i = 0; i <= n; ++i) {
        if (max < abs(d4_f(a + h * i))) max = abs(d4_f(a + h * i));
    }
    return max;
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
    for (int i = 1; i <= n; i++)
        integral = integral + 2.0 / 6.0 * h * f(a + h * i);
    return integral;
}

//////////////////////////////
// Оценка остаточных членов //
//////////////////////////////

// Теоретическая погрешность метода прямоугольников
double midpointAccuracy(double a, double b, double h) {
    return pow(h, 2) * (b - a) * max_d2_f(a, b, h) / 12.0;
}

// Теоретическая погрешность метода трапеций
double trapezoidalAccuracy(double a, double b, double h) {
    return pow(h, 2) * (b - a) * max_d2_f(a, b, h) / 6.0;
    // где max_d2_f это максимальная величина второй производной от f
}

// Теоретическая погрешность метода трех восьмых
double threeEighthsAccuracy(double a, double b, double h) {
    return pow(h, 4) * (b - a) * max_d4_f(a, b, h) / 80.0;
    // где max_d4_f это максимальная величина четвертой производной от f
}

// Теоретическая погрешность метода Симпсона
double simpsonAccuracy(double a, double b, double h) {
    return pow(h, 4) * (b - a) * max_d4_f(a, b, h) / 180.0;
    // где max_d4_f это максимальная величина четвертой производной от f
}

int main() {
    double a = 0.0;
    double b = M_PI;
    double h = 0.001;

    // Настройки точности вывода
    const auto digits = 10;
    cout << setw(digits);
    cout << fixed << std::setprecision(digits);

    // Результаты вычислений
    cout << "Midpoint sum " << midpointSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - midpointSum(a, b, h))
         << " meanwhile theoretical accuracy " << midpointAccuracy(a, b, h) << endl;
    cout << "Trapezoidal sum " << trapezoidalSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - trapezoidalSum(a, b, h))
         << " meanwhile theoretical accuracy " << trapezoidalAccuracy(a, b, h) << endl;
    cout << "Three eighths sum " << threeEighthSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - threeEighthSum(a, b, h))
         << " meanwhile theoretical accuracy " << threeEighthsAccuracy(a, b, h) << endl;
    cout << "Simpson sum " << simpsonSum(a, b, h) << " different from analytic "
         << fabs(analyticIntegral(a, b) - simpsonSum(a, b, h))
         << " meanwhile theoretical accuracy " << simpsonAccuracy(a, b, h) << endl;
    cout << "Analytic result " << analyticIntegral(a, b) << endl;
    return 0;
}
