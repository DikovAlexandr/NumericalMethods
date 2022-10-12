#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

double f_1(double x, double r, double o) {
    return 4.0 / 3.0 * 3.1415926535 * (x * x) * exp(-x) / r +
           4.0 / 15.0 * 3.1415926535 * (x * x * x * x) * exp(-x) / (r * r * r) * (3 * cos(o) * cos(o) - 1);
}

double f_2(double x, double r, double o) {
    return 4.0 / 3.0 * 3.1415926535 * (x + r) * exp(-r) +
           4.0 / 15.0 * 3.1415926535 * exp(-r) * r * r / (x + r) * (3 * cos(o) * cos(o) - 1);
}

double simpson(double a, double b, double r, double o) {
    double xi, xii, h, d, output;
    int n;
    n = 1001;
    h = ((b - a) / (n - 1));
    output = f_1(a, r, o) - f_1(b, r, o);
    for (int i = 1; i < n - 1; i + 2) {
        xi = a + i * h;
        xii = a + (i + 1) * h;
        output = output + (4.0 * f_1(xi, r, o)) + (2.0 * f_1(xii, r, o));
    }
    return (output * h) / 3;
}

double potential(double r, double o) {
    double I1, I2, I3, x1, x2, x3, x4, x5, x6, x7, w1, w2, w3, w4, w5, w6, w7, fi;
    x1 = 0.193044;
    x2 = 1.026665;
    x3 = 2.567877;
    x4 = 4.900353;
    x5 = 8.182153;
    x6 = 12.734180;
    x7 = 19.395728;

    w1 = 0.4093190;
    w2 = 0.4218313;
    w3 = 0.1471263;
    w4 = 0.02063351;
    w5 = 0.001074010;
    w6 = 0.00001586546;
    w7 = 0.00000003170315;
    I1 = simpson(0.0, r, r, o);
    I2 = w1 * f_2(x1, r, o) + w2 * f_2(x2, r, o) + w3 * f_2(x3, r, o) + w4 * f_2(x4, r, o) + w5 * f_2(x5, r, o) +
         w6 * f_2(x6, r, o) + w7 * f_2(x7, r, o);
    cout << "I1: " << I1 << endl;
    cout << "I2: " << I2 << endl;
    fi = (I1 + I2) * 3.0 / (8.0 * 3.1415926535);
    return fi;
}

int main() {

    return 0;
}

void graph() {
    double x0, xn, x, y, xi, xj, yi, c;
    int i, j, n, m = 1000;
    int parameter[m];
    double mass_1[m + 1], mass_2[m + 1], mass_3[m + 1], mass_4[m + 1], mass_5[m + 1], mass_6[m + 1];
    x0 = 0.01;
    xn = 100.0;
    ofstream fileOut;
    fileOut.open("data.txt");
    for (int k = 1; k < m; k++) {
        x = x0 + (k - 1) * (xn - x0) / (m - 1);
        mass_1[k] = x;
        mass_2[k] = potential(x, 0.0);
        mass_3[k] = potential(x, 3.1415926535 / 6);
        mass_4[k] = potential(x, 3.1415926535 / 4);
        mass_5[k] = potential(x, 3.1415926535 / 3);
        mass_6[k] = potential(x, 3.1415926535 / 2);
        fileOut << mass_1[m] << " " << mass_2[m] << " " << mass_3[m] << " " << mass_4[m] << " " << mass_5[m] << " "
                << mass_6[m];
    }
    fileOut.close();
}
