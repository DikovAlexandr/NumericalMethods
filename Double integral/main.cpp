
#include <iostream>
#include <cmath>

using namespace std;

///@brief define A and B semi-axes, inner circle radius R and angle of rotate PHI
#define A 1
#define B 0.5
#define R 1
#define PHI M_PI/4

double f(double x, double y) {
    #ifdef PHI
    double x1 = x * cos(PHI) - y * sin(PHI);
    double y1 = y * cos(PHI) + x * sin(PHI);
    x = x1; y = y1;
    #endif
    return pow(x, 2) * pow(y, 2);
}

double outerLimits(double x) {
    return sqrt(pow(R, 2) - pow(x, 2));
}

double innerLimits(double x) {
    #ifdef A
    if (fabs(x) > A)
        return 0;
    else return sqrt(1 - pow(x / A, 2)) * B;
    #endif
    #ifndef A
    return 0;
    #endif
}
/// @brief Gaussian 6-point integral, used in QC6 for double integral
/// @param c lower limit
/// @param d upper limit
/// @param x fixed x[i]
/// @return integral value
double QC6(double c, double d, double x) {
    double sum = 0;
    double mid = (c + d) / 2;
    double half = (d - c) / 2;
    double w[6] = {0.467913934572691, 0.467913934572691, 0.360761573048139, 0.360761573048139, 0.171324492379170, 0.171324492379170};
    double y[6] = {0.238619186083197, -0.238619186083197, 0.661209386466265, -0.661209386466265, 0.932469514203152, -0.932469514203152};
    for (int i = 0; i < 6; i++) 
        sum += half * w[i] * f(x, half * y[i] + mid);
    return sum;
}
/// @brief Gaussian 6-point double integral
/// @param a lower limit
/// @param b upper limit
/// @return double integral value
double QC6(double a, double b) {
    double sum = 0;
    double mid = (a + b) / 2;
    double half = (b - a) / 2;
    double w[6] = {0.467913934572691, 0.467913934572691, 0.360761573048139, 0.360761573048139, 0.171324492379170, 0.171324492379170};
    double x[6] = {0.238619186083197, -0.238619186083197, 0.661209386466265, -0.661209386466265, 0.932469514203152, -0.932469514203152};
    for (int i = 0; i < 6; i++)
        sum += half * w[i] * (QC6(innerLimits(half * x[i] + mid), outerLimits(half * x[i] + mid), half * x[i] + mid) + QC6(-outerLimits(half * x[i] + mid), -innerLimits(half * x[i] + mid), half * x[i] + mid));
    return sum;
}

/// @brief Gaussian 8-point integral, used in QC8 for double integral
/// @param c lower limit
/// @param d upper limit
/// @param x fixed x[i]
/// @return integral value
double QC8(double c, double d, double x) {
    double sum = 0;
    double mid = (c + d) / 2;
    double half = (d - c) / 2;
    double w[8] = {0.362683783378362, 0.362683783378362, 0.313706645877887, 0.313706645877887, 0.222381034453374, 0.222381034453374, 0.101228536290376, 0.101228536290376};
    double y[8] = {0.183434642495650, -0.183434642495650, 0.525532409916329, -0.525532409916329, 0.796666477413627, -0.796666477413627, 0.960289856497536, -0.960289856497536};
    for (int i = 0; i < 8; i++) 
        sum += half * w[i] * f(x, half * y[i] + mid);
    return sum;
}
/// @brief Gaussian 8-point double integral
/// @param a lower limit
/// @param b upper limit
/// @return double integral value
double QC8(double a, double b) {
    double sum = 0;
    double mid = (a + b) / 2;
    double half = (b - a) / 2;
    double w[8] = {0.362683783378362, 0.362683783378362, 0.313706645877887, 0.313706645877887, 0.222381034453374, 0.222381034453374, 0.101228536290376, 0.101228536290376};
    double x[8] = {0.183434642495650, -0.183434642495650, 0.525532409916329, -0.525532409916329, 0.796666477413627, -0.796666477413627, 0.960289856497536, -0.960289856497536};
    for (int i = 0; i < 8; i++)
        sum += half * w[i] * (QC8(innerLimits(half * x[i] + mid), outerLimits(half * x[i] + mid), half * x[i] + mid) + QC8(-outerLimits(half * x[i] + mid), -innerLimits(half * x[i] + mid), half * x[i] + mid));
    return sum;
}

/// @brief Gaussian 10-point integral, used in QC10 for double integral
/// @param c lower limit
/// @param d upper limit
/// @param x fixed x[i]
/// @return integral value
double QC10(double c, double d, double x) {
    double sum = 0;
    double mid = (c + d) / 2;
    double half = (d - c) / 2;
    double w[10] = {0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996, 0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688};
    double y[10] = {0.97390652851717174, 0.86506336668898454, 0.67940956829902444, 0.43339539412924721, 0.14887433898163122, -0.14887433898163119, -0.43339539412924721, -0.67940956829902444, -0.86506336668898454, -0.97390652851717174};
    for (int i = 0; i < 10; i++) 
        sum += half * w[i] * f(x, half * y[i] + mid);
    return sum;
}
/// @brief Gaussian 10-point double integral
/// @param a lower limit
/// @param b upper limit
/// @return double integral value
double QC10(double a, double b) {
    double sum = 0;
    double mid = (a + b) / 2;
    double half = (b - a) / 2;
    double w[10] = {0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996, 0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688};
    double x[10] = {0.97390652851717174, 0.86506336668898454, 0.67940956829902444, 0.43339539412924721, 0.14887433898163122, -0.14887433898163119, -0.43339539412924721, -0.67940956829902444, -0.86506336668898454, -0.97390652851717174};
    for (int i = 0; i < 10; i++)
        sum += half * w[i] * (QC10(innerLimits(half * x[i] + mid), outerLimits(half * x[i] + mid), half * x[i] + mid) + QC10(-outerLimits(half * x[i] + mid), -innerLimits(half * x[i] + mid), half * x[i] + mid));
    return sum;
}

int main() {
    #ifdef PHI
    cout << "Analytical answer: " << 85 * M_PI / 3072 << endl;
    #endif
	#ifndef PHI
	if (A && B)
		cout << "Analytical answer: " << 7 * M_PI / 192 << endl;
	else cout << "Analytical answer: " << M_PI / 24 << endl;;
	#endif
    cout << "6-point value: " << QC6(-R, R) << endl;
    cout << "8-point value: " << QC8(-R, R) << endl;
    cout << "10-point value: " << QC10(-R, R) << endl;
    return 0;
}
