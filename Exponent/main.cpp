#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <conio.h>
#include <limits>

using namespace std;

// Вычисление обратного факториала
double factorial(int k) {
    double out = 1;
    for (int i = 1; i <= k; ++i) {
        out /= i;
    }
    return out;
}

// Нахождение очередного слогаемого суммы
double fExp(double x, int k) {
    return pow(x, k) * factorial(k);
}

double func(int k) {
    return 1.0 / (pow(k, 4) * (pow(k, 2) + 1));
}

// Реализация метода Кэхэна для нахождения значения экспоненты
double kahanExp(double eps, double x, int &k) {
    double current;
    double element;
    double error = 0.0;
    double save = -1.0;
    double sum = 0.0;
    while (fabs(sum - save) > eps) {
        save = sum;
        element = fExp(x, k) - error;
        current = sum + element;
        error = current - sum - element;
        sum = current;
        k++;
    }
    return sum;
}

// Реализация метода Кэхэна для метода Куммера
double kahanKummer(double eps, const double standard, int &k) {
    double current;
    double element = 1.0;
    double error = 0.0;
    double sum = 0.0;
    k = 1;
    while (fabs(sum - standard) >= eps) {
        element = func(k) - error;
        current = sum + element;
        error = current - sum - element;
        sum = current;
        k++;
    }
    return sum;
}

// Тело функции
int main() {
    double eps = 1.0e-8;
    double x;
    int k = 0;
    cout << "Enter x" << endl;
    cin >> x;
    // Настройки точности вывода
    const auto digits = 8;
    cout << setw(digits);
    cout << fixed << setprecision(digits);
    // Нахождение значения экспоненты в опроделенной степени
    cout << "Summation result = " << kahanExp(eps, x, k) << endl;
    cout << "Analytics result = " << exp(x) << endl;
    cout << "N = " << k << endl;
    cout << "-------------------------------" << endl;
    // Нахождение суммы ряда методом Куммера
    const double answerAnalytic = 1.07667404746858;
    const double standard = answerAnalytic - pow(M_PI, 2) / 6.0 + pow(M_PI, 4) / 90.0;
    cout << "1/(k^2+1) series analysis" << endl;
    cout << "Summation result = " << pow(M_PI, 2) / 6.0 - pow(M_PI, 4) / 90.0 + kahanKummer(eps, standard, k)
         << endl;
    cout << "Analytics result = " << answerAnalytic << endl;
    cout << "N = " << k << endl;
}
