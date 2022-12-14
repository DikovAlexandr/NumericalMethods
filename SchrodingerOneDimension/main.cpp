#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <tuple>

using namespace std;

double f(double x, double y1, double y2) {
    return y2;
}

// Линейное слогаемное в уравнении Шредингера
double g(double x, double psi1, double psi2, double eigenEnergy) {
    return psi1 * 2.0 * (eigenEnergy - 10.0 / cosh(x));
}

//double g(double x, double psi1, double psi2, double eigenEnergy) {
//    return psi1 * 2.0 * (eigenEnergy - 10.0 / pow(cosh(x), 2));
//}

// Метод Кунге-Кутта для решения системы уравнений
pair<double, double> RK(double x, double psiN, double dPsiN, double h, double eigenEnergy) {
    double k1, k2, k3, k4, m1, m2, m3, m4;
    k1 = f(x, psiN, dPsiN);
    m1 = g(x, psiN, dPsiN, eigenEnergy);

    k2 = f(x + h / 2.0, psiN + h / 2.0 * k1, dPsiN + h / 2.0 * m1);
    m2 = g(x + h / 2.0, psiN + h / 2.0 * k1, dPsiN + h / 2.0 * m1, eigenEnergy);

    k3 = f(x + h / 2.0, psiN + h / 2.0 * k2, dPsiN + h / 2.0 * m2);
    m3 = g(x + h / 2.0, psiN + h / 2.0 * k2, dPsiN + h / 2.0 * m2, eigenEnergy);

    k4 = f(x + h, psiN + h * k3, dPsiN + h * m3);
    m4 = g(x + h, psiN + h * k3, dPsiN + h * m3, eigenEnergy);

    pair<double, double> ret;
    ret = make_pair(psiN + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4), dPsiN + h / 6.0 * (m1 + 2 * m2 + 2 * m3 + m4));
    return ret;
}

// Тело функции
int main() {
    const int N = 10000, M = 10000;
    double x0 = 0.005;
    double a, b;
    double hXLeft, hXRight, hEnergy;
    double psiN[3], dPsiN[3];
    // Началное и конечное значения энергии
    double aEnergy = 0.0, bEnergy = 10.0;
    double xn, xm;
    // Границы по x

    double xRight = 50, xLeft = -1 * xRight;
    double psiM[3], dPsiM[3];
    double E, W[M + 1];

    float psiLeft[N + 2][10], dPsiLeft[N + 2][10], psiRight[N + 2][10], dPsiRight[N + 2][10];
    double Energy[N + 2];
    double S1[N + 2], S2[N + 2];

    // Создаем файл для хранения значений энергии и вронскиана
    ofstream fout1;
    fout1.open(R"(..\\..\\SchrodingerOneDimension\\EnergyAndWronskian.csv)");

    // Нахождение собственных значений энергии
    hEnergy = (bEnergy - aEnergy) / M;
    for (int i = 0; i <= M; ++i) {
        ////////////////////
        // Левая половина //
        ////////////////////

        b = x0;
        a = xLeft;
        hXLeft = (b - a) / N;
        E = aEnergy + hEnergy * i;

        // Считаем что волновая функция имеет асимптотический вид
        psiN[2] = exp(xLeft * sqrt(2.0 * E));
        // Производная от волновой функции
        dPsiN[2] = sqrt(2.0 * E) * exp(xLeft * sqrt(2.0 * E));

        // Ведем асимптотическое решение к нулю
        xn = xLeft;
        for (int k = 0; k <= N; ++k) {
            psiN[1] = psiN[2];
            dPsiN[1] = dPsiN[2];
            psiN[2] = RK(xn, psiN[1], dPsiN[1], hXLeft, E).first;
            dPsiN[2] = RK(xn, psiN[1], dPsiN[1], hXLeft, E).second;
            xn = xLeft + hXLeft * k;
        }

        /////////////////////
        // Правая половина //
        /////////////////////

        b = xRight;
        a = x0;
        hXLeft = (a - b) / N;

        psiM[2] = exp(-xRight * sqrt(2.0 * E));
        dPsiM[2] = -sqrt(2.0 * E) * exp(-xRight * sqrt(2.0 * E));

        // Ведем асимптотическое решение к нулю
        xm = xRight;
        for (int k = 0; k <= N; ++k) {
            psiM[1] = psiM[2];
            dPsiM[1] = dPsiM[2];
            psiM[2] = RK(xm, psiM[1], dPsiM[1], hXLeft, E).first;
            dPsiM[2] = RK(xm, psiM[1], dPsiM[1], hXLeft, E).second;
            xm = xRight + hXLeft * k;
        }

        // Если значение энергии является особственным функции будут сшиваться
        // Находим определитель Вронского
        W[i] = dPsiN[2] * psiM[2] - psiN[2] * dPsiM[2];
        fout1 << E << " " << W[i] << endl;
    }
    fout1.close();

    // Отбор энергий
    double energyI = 0.0;
    int n = 0;
    for (int i = 1; i <= M; ++i) {
        if (W[i] * W[i - 1] <= 0) {
            energyI = energyI + hEnergy;
            Energy[n] = energyI;
            cout << Energy[n] << endl;
            n++;
        } else energyI = energyI + hEnergy;
    }

    for (int i = 0; i < n; ++i) {
        ///////////////////
        // Функция слева //
        ///////////////////

        E = Energy[i];
        xRight = -xLeft;

        b = x0;
        a = xLeft;
        hXLeft = (b - a) / N;

        // Волновая функция и ее производная на асимптотиках
        psiLeft[0][i] = exp(xLeft * sqrt(2.0 * Energy[i]));
        dPsiLeft[0][i] = sqrt(2.0 * Energy[n]) * exp(xLeft * sqrt(2.0 * Energy[i]));
        xn = xLeft;

        for (int k = 0; k < N; ++k) {
            psiLeft[k + 1][i] = RK(xn, psiLeft[k][i], dPsiLeft[k][i], hXLeft, Energy[i]).first;
            dPsiLeft[k + 1][i] = RK(xn, psiLeft[k][i], dPsiLeft[k][i], hXLeft, Energy[i]).second;
            xn = xLeft + hXLeft * k;
        }

        psiLeft[N][i] = RK(xn, psiLeft[N - 1][i], dPsiLeft[N - 1][i], hXLeft, Energy[i]).first;

        // Нормировка методом Симпсона
        S1[i] = pow(psiLeft[0][i], 2) - pow(psiLeft[N][i], 2);
        for (int j = 1; j < N; j += 2) {
            S1[i] = S1[i] + 4.0 * pow(psiLeft[j][i], 2) + 2.0 * pow(psiLeft[j + 1][i], 2);
        }
        S1[i] = S1[i] * hXLeft / 3.0;

        ////////////////////
        // Функция справа //
        ////////////////////

        b = xRight;
        a = x0;
        hXRight = (a - b) / N;

        // Волновая функция и ее производная на асимптотиках
        psiRight[0][i] = exp(-xRight * sqrt(2.0 * Energy[i]));
        dPsiRight[0][i] = -sqrt(2.0 * Energy[n]) * exp(-xRight * sqrt(2.0 * Energy[i]));
        xm = xRight;

        for (int k = 0; k < N; ++k) {
            psiRight[k + 1][i] = RK(xm, psiRight[k][i], dPsiRight[k][i], hXRight, Energy[i]).first;
            dPsiRight[k + 1][i] = RK(xm, psiRight[k][i], dPsiRight[k][i], hXRight, Energy[i]).second;
            xm = xRight + hXRight * k;
        }

        psiRight[N][i] = RK(xm, psiRight[N - 1][i], dPsiRight[N - 1][i], hXRight, Energy[i]).first;

        // Нормировка методом Симпсона
        S2[i] = -pow(psiRight[0][i], 2) + pow(psiRight[N][i], 2);
        for (int j = 1; j < N; ++j) {
            S2[i] = S2[i] + 4.0 * pow(psiRight[j][i], 2) + 2.0 * pow(psiRight[j + 1][i], 2);
        }
        S2[i] = -1 * S2[i] * hXRight / 3.0;
    }

    // Значения варьированной волновой функции слева
    ofstream fout2;
    fout2.open(R"(..\\..\\SchrodingerOneDimension\\Psi.csv)");
    int sign;
    for (int i = 0; i <= N; ++i) {
        xn = xLeft + hXLeft * i;
        fout2 << xn;
        for (int j = 0; j < n; ++j) {
            //psiLeft[i][j] = pow(psiLeft[i][j], 2) / (S1[j] + S2[j]);
            psiLeft[i][j] = psiLeft[i][j] / sqrt(S1[j] + S2[j]);
            if ((n - j + 1)%2 != 0) sign = -1;
            else sign = 1;
            fout2 << " " << sign * psiLeft[i][j];
            //fout2 << " " << sign * dPsiLeft[i][j];
        }
        fout2 << endl;
    }

    // Значения варьированной волновой функции справа
    for (int i = 0; i <= N; ++i) {
        xm = x0 + (-1) * hXRight * i;
        fout2 << xm;
        for (int j = 0; j < n; ++j) {
            // psiRight[i][j] = pow(psiRight[i][j], 2) / (S1[j] + S2[j]);
            psiRight[N - i][j] = psiRight[N - i][j] / sqrt(S1[j] + S2[j]);
            fout2 << " " << psiRight[N - i][j];
            //fout2 << " " << dPsiRight[N - i][j];
        }
        fout2 << endl;
    }
    fout2.close();
}