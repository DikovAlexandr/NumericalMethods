#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Предварительно получанное разложение потенциала по шаровым функциям разбиваем на два слогаемых
// Первая функция содержит подытегральное выражение от слогаемого с собственным интегралом
double f_1(double x, double r, double o) {
    return 4.0 / 3.0 * M_PI * pow(x, 2) * exp(-x) / r +
           4.0 / 15.0 * M_PI * pow(x, 4) * exp(-x) / pow(r, 3) * (3 * pow(cos(o), 2) - 1);
}

// Вторая функция содержит подытегральное выражение от слогаемого несобственным интегралом
double f_2(double x, double r, double o) {
    return 4.0 / 3.0 * M_PI * (x + r) * exp(-r) +
           4.0 / 15.0 * M_PI * exp(-r) * pow(r, 2) / (x + r) * (3 * pow(cos(o), 2) - 1);
}

// Метод Симпсона для собственного интеграла
double simpson(double a, double b, double r, double o) {
    double xi, xii, h, d, output;
    int n = 1000;
    h = ((b - a) / (n - 1));
    output = f_1(a, r, o) - f_1(b, r, o);
    for (int i = 1; i < n; i += 2) {
        xi = a + i * h;
        xii = a + (i + 1) * h;
        output = output + (4.0 * f_1(xi, r, o)) + (2.0 * f_1(xii, r, o));
    }
    return (output * h) / 3;
}

// Вызов Симпсона и метод Гаусса-Кристоффеля с весовой функцией exp(-x) для несобсвенного интеграла
double potential(double r, double o) {
    double I1, I2, fi;
    double x[7] = {0.193044, 1.026665, 2.567877, 4.900353, 8.182153, 12.734180, 19.395728};
    double w[7] = {0.4093190, 0.4218313, 0.147126, 0.02063351, 0.001074010, 0.00001586546, 0.00000003170315};
    I1 = simpson(0.0, r, r, o);
    I2 = w[0] * f_2(x[0], r, o) + w[1] * f_2(x[1], r, o) + w[2] * f_2(x[2], r, o) + w[3] * f_2(x[3], r, o) +
         w[4] * f_2(x[4], r, o) + w[5] * f_2(x[5], r, o) + w[6] * f_2(x[6], r, o);
    fi = (I1 + I2) * 3.0 / (8.0 * M_PI);
    return fi;
}

// Тело функции
int main() {
    double x0, xn, x;
    int m = 1000;
    double xArray[m], theta0[m], theta30[m], theta45[m], theta60[m], theta90[m];
    x0 = 0.01;
    xn = 100.0;
    ofstream fout;
    fout.open(R"(..\\..\\NumericalMethods\\Potencial of field\\potential.csv)");
    for (int k = 1; k < m; k++) {
        x = x0 + (k - 1) * (xn - x0) / (m - 1);
        xArray[k] = x;
        theta0[k] = potential(x, 0.0);
        theta30[k] = potential(x, 3.1415926535 / 6);
        theta45[k] = potential(x, 3.1415926535 / 4);
        theta60[k] = potential(x, 3.1415926535 / 3);
        theta90[k] = potential(x, 3.1415926535 / 2);
        fout << xArray[k] << " " << theta0[k] << " " << theta30[k] << " " << theta45[k] << " " << theta60[k] << " "
             << theta90[k] << endl;
    }
    fout.close();
    return 0;
}