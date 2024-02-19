#include <iostream>
#include <limits>
#include <cmath>

using namespace std;

// Машинный эпсилон double
double getEpsilonDouble() {
    double eps = 0.5;
    while (1.0 + eps != 1.0) {
        eps *= 0.5;
    }
    return eps;
};

// Машинный эпсилон float
float getEpsilonFloat() {

    float eps = 0.5;
    while (1.0 + eps != 1.0) {
        eps *= 0.5;
    }
    return eps;
};

// Машинный ноль double
double getMachineZeroDouble() {
    double zero = 1.0, save;
    do {
        save = zero;
        zero /= 2;
    } while (zero / 2 != 0);
    return save;
};

// Машинный ноль float
float getMachineZeroFloat() {
    float zero = 1.0, save;
    do {
        save = zero;
        zero /= 2;
    } while (zero / 2 != 0);
    return save;
};

// Бесконечность double
double getInfinityDouble() {
    double inf = 1.0, k = 2.0;
    while (!isinf(inf * k)) {
        inf *= k;
        if (isinf(inf * k)) k = 1 + (k - 1) / 2;
        if (fabs(1 - k) < 0.01) break;
    }
    return inf;
};

// Бесконечность float
float getInfinityFloat() {
    float inf = 1.0, k = 2.0;
    while (!isinf(inf * k)) {
        inf *= k;
        if (isinf(inf * k)) k = 1 + (k - 1) / 2;
        if (fabs(1 - k) < 0.01) break;
    }
    return inf;
};

int main() {
    cout << "Machine zero double " << getMachineZeroDouble() << " // " << numeric_limits<double>::denorm_min() << endl
         << "Machine zero float " << getMachineZeroFloat() << " // " << numeric_limits<float>::denorm_min() << endl
         << "Epsilon double " << getEpsilonDouble() << " // " << numeric_limits<double>::epsilon() << endl
         << "Epsilon float " << getEpsilonFloat() << " // " << numeric_limits<float>::epsilon() << endl
         << "Infinity double " << getInfinityDouble() << " // " << numeric_limits<double>::max() << endl
         << "Infinity float " << getInfinityFloat() << " // " << numeric_limits<float>::max() << endl;
    return 0;
}
