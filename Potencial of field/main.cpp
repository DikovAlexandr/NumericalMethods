#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double f_1(double x, double r, double o) {
    return 4.0 / 3.0 * M_PI * pow(x, 2) * exp(-x) / r +
           4.0 / 15.0 * M_PI * pow(x, 4) * exp(-x) / pow(r, 3) * (3 * pow(cos(o), 2) - 1);
}

double f_2(double x, double r, double o) {
    return 4.0 / 3.0 * M_PI * (x + r) * exp(-r) +
           4.0 / 15.0 * M_PI * exp(-r) * pow(r, 2) / (x + r) * (3 * pow(cos(o), 2) - 1);
}

double simpson(double a, double b, double r, double o) {
    double xi, xii, h, d, output;
    int n = 1000;
    h = ((b - a) / (n - 1));
    output = f_1(a, r, o) - f_1(b, r, o);
    for (int i = 1; i < n ; i += 2) {
        xi = a + i * h;
        xii = a + (i + 1) * h;
        output = output + (4.0 * f_1(xi, r, o)) + (2.0 * f_1(xii, r, o));
    }
    return (output * h) / 3;
}

double potential(double r, double o) {
    double I1, I2, I3, x2, x3, x4, x5, x6, x7, w1, w2, w3, w4, w5, w6, w7, fi;
    double x[7] = {0.193044, 1.026665, 2.567877, 4.900353, 8.182153, 12.734180, 19.395728};
    double w[7] = {0.4093190, 0.4218313, 0.147126, 0.02063351, 0.001074010, 0.00001586546, 0.00000003170315};
    I1 = simpson(0.0, r, r, o);
    I2 = w[0] * f_2(x[0], r, o) + w[1] * f_2(x[1], r, o) + w[2] * f_2(x[2], r, o) + w[3] * f_2(x[3], r, o) + w[4] * f_2(x[4], r, o) +
         w[5] * f_2(x[5], r, o) + w[6] * f_2(x[6], r, o);
    fi = (I1 + I2) * 3.0 / (8.0 * M_PI);
    return fi;
}

int main() {
    double x0, xn, x, y, xi, xj, yi, c;
    int i, j, n, m = 1000;
    double mass_1[m + 1], mass_2[m + 1], mass_3[m + 1], mass_4[m + 1], mass_5[m + 1], mass_6[m + 1];
    x0 = 0.01;
    xn = 100.0;
    ofstream fout;
    fout.open(R"(..\\..\\NumericalMethods\\Potencial of field\\potential.csv)");
    for (int k = 1; k < m; k++) {
        x = x0 + (k - 1) * (xn - x0) / (m - 1);
        mass_1[k] = x;
        mass_2[k] = potential(x, 0.0);
        mass_3[k] = potential(x, 3.1415926535 / 6);
        mass_4[k] = potential(x, 3.1415926535 / 4);
        mass_5[k] = potential(x, 3.1415926535 / 3);
        mass_6[k] = potential(x, 3.1415926535 / 2);
        fout << mass_1[k] << " " << mass_2[k] << " " << mass_3[k] << " " << mass_4[k] << " " << mass_5[k] << " "
                << mass_6[k] << endl;
    }
    fout.close();
    return 0;
}