#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// Разобъем интеграл на два
// sin(x)/x интегрируется от -2 до 3
double f1(double x) {
    return sin(x) / x;
}

// cos(x)/x интегрирутся от 2 до 3
double f2(double x) {
    return cos(x) / x;
}

// 10-точечный метод Гаусса
double QC10(double c, double d, double f(double x)) {
    double sum = 0;
    double mid = (c + d) / 2;
    double half = (d - c) / 2;
    double w[10] = {0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996, 0.295524224714753,
                    0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688};
    double y[10] = {0.97390652851717174, 0.86506336668898454, 0.67940956829902444, 0.43339539412924721,
                    0.14887433898163122, -0.14887433898163119, -0.43339539412924721, -0.67940956829902444,
                    -0.86506336668898454, -0.97390652851717174};
    for (int i = 0; i < 10; i++)
        sum += half * w[i] * f(half * y[i] + mid);
    return sum;
}

// Метод Симпсона
double simpsonSumNew(double a, double b, double n, double f(double x)) {
    double h = (b - a) / (n * 2.0);
    double integral = 0, point;
    for (int i = 0; i < n; ++i) {
        point = a + i * h * 2.0;
        integral += f(point) + 4.0 * f(point + h) + f(point + 2.0 * h);
    }
    return integral * h / 3.0;
}

// Тело функции
int main() {
    // Настройки точности вывода
    const auto digits = 5;
    cout << setw(digits);
    cout << fixed << std::setprecision(digits);
    // Опеределяем константы
    int n = 1000;
    double I1, I2, S3;
    double a = -2.0, b = 3.0;
    // Вычисляем интегралы
    I1 = QC10(a, b, f1);
    a = 2.0;
    I2 = simpsonSumNew(a, b, n, f2);
    cout << "Numerical result " << I1 + I2 << ", analytic result " << 3.45407 - 0.30335;
}
