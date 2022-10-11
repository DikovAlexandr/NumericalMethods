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
double getEpsilonFloat() {

    float eps = 0.5;
    while (1.0 + eps != 1.0) {
        eps *= 0.5;
    }
    return eps;
};
// Машинный ноль double
double getMachineZeroDouble() {
    double zero = 1.0, save;
    while (zero * 0.5 != 0) {
        save = zero;
        zero *= 0.5;
    }
    return save;
};
// Машинный ноль float
double getMachineZeroFloat() {
    float zero = 1.0, save;
    while (zero * 0.5 != 0) {
        save = zero;
        zero *= 0.5;
    }
    return save;
};
// Бесконечность double
double getInfinityDouble() {
    double inf = 1.0;
    while (!isinf(inf * 2)) {
        inf *= 2;
    }
    return inf;
};
// Бесконечность float
double getInfinityFloat() {
    float inf = 1.0;
    while (!isinf(inf * 2)) {
        inf *= 2;
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
