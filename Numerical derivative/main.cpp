#include <iostream>
#include <fstream>
#include <cmath>
// программа которая вычисляет производную первого и второго порядка численно и аналитически от произвольной функции
using namespace std;
// Исходная функция
double f(double x) {
    return exp(-x) * pow(x, 2);
}
// Аналитическая первая производная
double f1(double x) {
    return (2 * x - pow(x, 2)) * exp(-x);
}
// Аналитическая вторая производная
double f2(double x) {
    return (2 - 4 * x + pow(x, 2)) * exp(-x);
}
// Численная первая производная
double diff1(double x, double h) {
    double fl, fr, fc;
    // Приближенно вычисляем первую производную различными способами
    //fl = (f(x) - f(x - h)) / h; // левая
    //fr = (f(x + h) - f(x)) / h; // правая
    fc = (f(x + h) - f(x - h)) / (2 * h); // центральная
    return fc;
}
// Численная вторая производная
double diff2(double x, double h) {
    double f2;
    // Приближенно вычисляем вторую производную
    f2 = (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
    return f2;
}
int main() {
    ofstream fout;
    fout.open("plotA.txt");
    double x, a, b, h, i;
    cout << "Enter range and step:\n";
    //cin >> a >> b >> h;
    a = 0;
    b = 1;
    h = 0.000001;
    i = a;
    while (i < b) {
        fout << i << " " << f(i) << " " << f1(i) << " " << diff1(i, h) << endl;
        i += h;
    }
    fout.close();
    i = a;
    fout.open("plotB.txt");
    while (i < b) {
        fout << i << " " << f(i) << " " << f2(i) << " " << diff2(i, h) << endl;
        i += h;
    }
    fout.close();
    return 0;
}
