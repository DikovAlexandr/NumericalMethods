#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

// ?
double f(double x, double y, double z) {
    return z;
}

// Первый потенциал U0/ch(x)
double g1(double x, double y, double z, double E) {
    double U0, U;
    //Задание потенциала
    U0 = 1.0;
    U = -1.0 * U0 / cosh(x);
    return -2.0 * (E - U) * y;
}

// Второй потенциал U0/cosh(x)**2
double g2(double x, double y, double z, double E) {
    double U0, U;
    //Задание потенциала
    U0 = 1.0;
    U = -1.0 * U0 / pow(cosh(x), 2);
    return -2.0 * (E - U) * y;
}

tuple<double, double, double> algorithm(double x, double y, double z, double E, double h, double (*g)(double, double, double, double)) {
    double k1, l1, k2, l2, k3, l3, k4, l4;
    k1 = f(x, y, z);
    l1 = g(x, y, z, E);
    k2 = f(x + h / 2.0, y + h / 2.0 * k1, z + h / 2.0 * l1);
    l2 = g(x + h / 2.0, y + h / 2.0 * k1, z + h / 2.0 * l1, E);
    k3 = f(x + h / 2.0, y + h / 2.0 * k2, z + h / 2.0 * l2);
    l3 = g(x + h / 2.0, y + h / 2.0 * k2, z + h / 2.0 * l2, E);
    k4 = f(x + h, y + h * k3, z + h * l3);
    l4 = g(x + h, y + h * k3, z + h * l3, E);
    y = y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    z = z + h / 6.0 * (l1 + 2.0 * l2 + 2.0 * l3 + l4);
    x = x + h;
    return make_tuple(x, y, z);
}

// Строим численное решение уравнения Шредингера методом Рунге-Кутты 4-го порядка
pair<double, double>
RK4(double x0, double y0, double z0, double h, double xs, double (*g)(double, double, double, double), double E) {
    double k1, k2, k3, k4;
    double l1, l2, l3, l4;
    double x, y, z;
    x = x0;
    y = y0;
    z = z0;
    // xs это точка сшивки
    if (h > 0.0) {
        while (true){
            if (x + h > xs) {
                k1 = f(x, y, z);
                l1 = g(x, y, z, E);
                k2 = f((x + xs) / 2.0, y + (xs - x) / 2.0 * k1, z + (xs - x) / 2.0 * l1);
                l2 = g((x + xs) / 2.0, y + (xs - x) / 2.0 * k1, z + (xs - x) / 2.0 * l1, E);
                k3 = f((x + xs) / 2.0, y + (xs - x) / 2.0 * k2, z + (xs - x) / 2.0 * l2);
                l3 = g((x + xs) / 2.0, y + (xs - x) / 2.0 * k2, z + (xs - x) / 2.0 * l3, E);
                k4 = f(xs, y + (xs - x) * k3, z + (xs - x) * l3);
                l4 = g(xs, y + (xs - x) * k3, z + (xs - x) * l3, E);
                y = y + (xs - x) / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                z = z + (xs - x) / 6.0 * (l1 + 2.0 * l2 + 2.0 * l3 + l4);
                return make_pair(y, z);
            } else if (x + h == xs) {
                tuple<double, double, double> temp = algorithm(x, y, z, E, h, g);
                return make_pair(get<1>(temp), get<2>(temp));
            } else {
                tuple<double, double, double> temp = algorithm(x, y, z, E, h, g);
                x = get<1>(temp);
                y = get<1>(temp);
                z = get<2>(temp);
            }
        }
    } else {
        while (true){
            if (x + h < xs) {
                k1 = f(x, y, z);
                l1 = g(x, y, z, E);
                k2 = f((x + xs) / 2.0, y + (xs - x) / 2.0 * k1, z + (xs - x) / 2.0 * l1);
                l2 = g((x + xs) / 2.0, y + (xs - x) / 2.0 * k1, z + (xs - x) / 2.0 * l1, E);
                k3 = f((x + xs) / 2.0, y + (xs - x) / 2.0 * k2, z + (xs - x) / 2.0 * l2);
                l3 = g((x + xs) / 2.0, y + (xs - x) / 2.0 * k2, z + (xs - x) / 2.0 * l3, E);
                k4 = f(xs, y + (xs - x) * k3, z + (xs - x) * l3);
                l4 = g(xs, y + (xs - x) * k3, z + (xs - x) * l3, E);
                y = y + (xs - x) / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
                z = z + (xs - x) / 6.0 * (l1 + 2.0 * l2 + 2.0 * l3 + l4);
                return make_pair(y, z);
            } else if (x + h == xs) {
                tuple<double, double, double> temp = algorithm(x, y, z, E, h, g);
                return make_pair(get<1>(temp), get<2>(temp));
            } else {
                tuple<double, double, double> temp = algorithm(x, y, z, E, h, g);
                x = get<1>(temp);
                y = get<1>(temp);
                z = get<2>(temp);
            }
        }
    }
}

// Сшивка функций
double stitch(double E, double (*g)(double, double, double, double)) {
    double xl, xr, xs;
    double kappa, h;
    double yl, zl, yr, zr;
    double yfl, zfl, yfr, zfr;
    // Начальные точки
    xl = -10.0;
    xr = 10.0;
    // Точка сшивки
    xs = 1e-3;
    // Шаг RK4
    h = 1e-5;
    kappa = sqrt(2.0 * fabs(E));
    // Начальные значения волновых функций и их производных слева и справа
    yl = exp(kappa * xl);
    yr = exp(-kappa * xr);
    zl = kappa * exp(kappa * xl);
    zr = -kappa * exp(-kappa * xr);
    pair<double, double> left = RK4(xl, yl, zl, h, xs, g, E);
    pair<double, double> right = RK4(xr, yr, zr, -h, xs, g, E);
    return left.second * right.first - right.second * left.first;
}

// Нахождение зависимости вранскиана от E
void graph(double E1, double E2, double (*g)(double, double, double , double)) {
    double h, E, det;
    // Шаг для зависимости Вронскиана от энергии
    h = 1e-3;
    E = E1;
    ofstream fout;
    fout.open(R"(..\\..\\NumericalMethods\\Schrodinger\\Wronski.csv)");
    while (E < E2) {
        fout << E << " " << stitch(E, g) << endl;
        E = E + h;
    }
    fout.close();
}

void energy(double (*g)(double, double, double, double)) {
    double E0, det0;
    double E1, det1;
    int star;
    star = 0;
    ofstream fout;
    ifstream in(R"(..\\..\\NumericalMethods\\Schrodinger\\Wronski.csv)");
    in >> E0 >> det0;
    while (!in.eof()) {
        in  >> E1 >> det1;
        if (det0 * det1 < 0) cout << (E0 + E1) / 2.0;
        else if (det0 == 0.0) {
            cout << E0;
        } else if (det1 == 0.0) {
            cout << E1;
        }
        E0 = E1;
        det0 = det1;
    }
}

// Тело функции
int main() {
    double E1, E2;
    E1 = -2.0;
    E2 = 0.0;
    graph(E1, E2, g2);
    energy(g2);
    return 0;
}