#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

double f(double x, double y1, double y2) {
    return y2;
}

double g(double x, double y1, double y2, double E_pot) {
    return y1 * 2.0 * (E_pot - 1.0 / cosh(x));
}

double RK1(double x, double y1, double y2, double h, double E_pot) {
    double m1, m2, m3, m4, k1, k2, k3, k4, der_value;
    k1 = f(x, y1, y2);
    m1 = g(x, y1, y2, E_pot);

    k2 = f(x + h / 2.0, y1 + h / 2.0 * k1, y2 + h / 2.0 * m1);
    m2 = g(x + h / 2.0, y1 + h / 2.0 * k1, y2 + h / 2.0 * m1, E_pot);

    k3 = f(x + h / 2.0, y1 + h / 2.0 * k2, y2 + h / 2.0 * m2);
    m3 = g(x + h / 2.0, y1 + h / 2.0 * k2, y2 + h / 2.0 * m2, E_pot);

    k4 = f(x + h, y1 + h * k3, y2 + h * m3);
    m4 = g(x + h, y1 + h * k3, y2 + h * m3, E_pot);

    return y1 + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
}

double RK2(double x, double y1, double y2, double h, double E_pot) {
    double m1, m2, m3, m4, k1, k2, k3, k4, der_value;

    k1 = f(x, y1, y2);
    m1 = g(x, y1, y2, E_pot);

    k2 = f(x + h / 2.0, y1 + h / 2.0 * k1, y2 + h / 2.0 * m1);
    m2 = g(x + h / 2.0, y1 + h / 2.0 * k1, y2 + h / 2.0 * m1, E_pot);

    k3 = f(x + h / 2.0, y1 + h / 2.0 * k2, y2 + h / 2.0 * m2);
    m3 = g(x + h / 2.0, y1 + h / 2.0 * k2, y2 + h / 2.0 * m2, E_pot);

    k4 = f(x + h, y1 + h * k3, y2 + h * m3);
    m4 = g(x + h, y1 + h * k3, y2 + h * m3, E_pot);

    return y2 + h / 6.0 * (m1 + 2 * m2 + 2 * m3 + m4);
}

int main() {
    int N = 1000, M = 1000, M1 = 10;
    double x0 = 1e-6, a, b, h1, h2, h3, y1n[3], xn, a_e = 0.0, b_e = 1.0;
    double y2n[3], x_left = -0.5e2, x_right, E, c, d, u = 1e-6, x, p;
    double y1m[3], y2m[3], xm, W[M + 1], y1[N + 2], y2[N + 1], y3[N + 2], y4[N + 1], y_h, Energy[2];
    double S1 = 0.0, S2 = 0.0;
    int k, i;

    h3 = (b_e - a_e) / M;
    ofstream fout1;
    fout1.open(R"(myfile1.txt)"); // 'myfile1.txt'
    for (i = 0; i <= M; ++i) {
        x_right = -x_left;
        b = x0;
        a = x_left;
        h1 = (b - a) / N;
        E = a_e + h3 * i;

        y1n[2] = exp(x_left * sqrt(2.0 * E));
        y2n[2] = sqrt(2.0 * E) * exp(x_left * sqrt(2.0 * E));
        xn = x_left;

        for (k = 1; k <= N; ++k) {
            y1n[1] = y1n[2];
            y2n[1] = y2n[2];
            y1n[2] = RK1(xn, y1n[1], y2n[1], h1, E);
            y2n[2] = RK2(xn, y1n[1], y2n[1], h1, E);
            xn = x_left + h1 * k;
        }

        b = x_right;
        a = x0;
        h2 = (a - b) / N;

        y1m[2] = exp(-x_right * sqrt(2.0 * E));
        y2m[2] = -sqrt(2.0 * E) * exp(-x_right * sqrt(2.0 * E));
        xm = x_right;

        for (k = 1; k <= N; ++k) {
            y1m[1] = y1m[2];
            y2m[1] = y2m[2];
            y1m[2] = RK1(xm, y1m[1], y2m[1], h2, E);
            y2m[2] = RK2(xm, y1m[1], y2m[1], h2, E);
            xm = x_right + h2 * k;
        }
        W[i] = y2n[2] * y1m[2] - y1n[2] * y2m[2];
        fout1 << -E << " " << W[i] << endl;
    }
    fout1.close();

    d = 0.0;
    k = 1;
    for (i = 2; i <= M; ++i) {
        if (W[i] * W[i - 1] <= 0) {
            d = d + h3;
            Energy[k] = d;
            k = k + 1;
        } else d = d + h3;
    }

    cout << Energy << endl;
    x_right = -x_left;
    b = x0;
    a = x_left;
    h1 = (b - a) / N;
    k = 1;
    E = Energy[k];
    y1[1] = exp(x_left * sqrt(2.0 * E));
    y2[1] = sqrt(2.0 * E) * exp(x_left * sqrt(2.0 * E));
    xn = x_left;

    for (k = 1; k < N; ++k) {
        y1[k + 1] = RK1(xn, y1[k], y2[k], h1, E);
        y2[k + 1] = RK2(xn, y1[k], y2[k], h1, E);
        xn = x_left + h1 * k;
    }

    y1[N + 1] = RK1(xn, y1[N], y2[N], h1, E);
    S1 = pow(y1[1], 2) - pow(y1[N + 1], 2);

    for (i = 2; i <= N; i += 2) {
        S1 = S1 + 4.0 * pow(y1[i], 2) + 2.0 * pow(y1[i + 1], 2);
    }

    S1 = S1 * h1 / 3.0;
    cout << S1 << endl;

    b = x_right;
    a = x0;
    h2 = (a - b) / N;
    y3[1] = exp(-x_right * sqrt(2.0 * E));
    y4[1] = -sqrt(2.0 * E) * exp(-x_right * sqrt(2.0 * E));
    xm = x_right;

    for (k = 1; k < N; ++k) {
        y3[k + 1] = RK1(xm, y3[k], y4[k], h2, E);
        y4[k + 1] = RK2(xm, y3[k], y4[k], h2, E);
        xm = x_right + h2 * k;
    }

    y3[N + 1] = RK1(xm, y3[N], y4[N], h2, E);
    S2 = -pow(y3[1], 2) + pow(y3[N + 1], 2);

    for (i = 2; i <= N; ++i) {
        S2 = S2 + 4.0 * pow(y3[i], 2) + 2.0 * pow(y3[i + 1], 2);
    }

    S2 = -S2 * h2 / 3.0;
    cout << S2 << endl;

    ofstream fout2;
    fout2.open(R"(..\\..\\NumericalMethods\\Schrodinger\\func11.csv)"); // "func11.txt"
    b = x0;
    a = x_left;
    for (i = 1; i <= N + 1; ++i) {
        y1[i] = pow(y1[i], 2) / (S1 + S2);
        xn = x_left + h1 * (i - 1);
        fout2 << xn << " " << y1[i] << endl;
    }
    fout2.close();

    ofstream fout3;
    fout3.open(R"(..\\..\\NumericalMethods\\Schrodinger\\func12.csv)"); // "func12.txt"
    b = x_right;
    a = x0;
    for (i = 1; i <= N + 1; ++i) {
        y3[i] = pow(y3[i], 2) / (S1 + S2);
        xm = x_right + h2 * (i - 1);
        fout3 << xm << " " << y3[i] << endl;
    }
    fout3.close();
}