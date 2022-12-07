#include <iostream>
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
    return psi1 * 2.0 * (eigenEnergy - 1.0 / cosh(x));
}

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
    int N = 1000, M = 1000, M1 = 10;
    double x0 = 1e-6;
    double a, b;
    double hX, h2, hEnergy;
    double psiN[3], dPsiN[3];
    // Началное и конечное значения энергии
    double aEnergy = 0.0, bEnergy = 1.0;
    double xn, xm;
    // Границы по x
    double xLeft = -0.5e2, xRight = -1 * xLeft;
    double psiM[3], dPsiM[3];
    double E, W[M + 1];
    double psiLeft[N + 2], dPsiLeft[N + 1], psiRight[N + 2], dPsiRight[N + 1];
    double Energy[2];
    double S1 = 0.0, S2 = 0.0;

    // Создаем файл для хранения значений энергии и вронскиана
    ofstream fout1;
    fout1.open(R"(EnergyAndWronskian.txt)");

    // Нахождение собственных значений энергии
    hEnergy = (bEnergy - aEnergy) / M;
    for (int i = 0; i <= M; ++i) {
        ////////////////////
        // Левая половина //
        ////////////////////

        b = x0;
        a = xLeft;
        hX = (b - a) / N;
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
            psiN[2] = RK(xn, psiN[1], dPsiN[1], hX, E).first;
            dPsiN[2] = RK(xn, psiN[1], dPsiN[1], hX, E).second;
            xn = xLeft + hX * k;
        }

        /////////////////////
        // Правая половина //
        /////////////////////

        b = xRight;
        a = x0;
        hX = (a - b) / N;

        psiM[2] = exp(-xRight * sqrt(2.0 * E));
        dPsiM[2] = -sqrt(2.0 * E) * exp(-xRight * sqrt(2.0 * E));

        // Ведем асимптотическое решение к нулю
        xm = xRight;
        for (int k = 0; k <= N; ++k) {
            psiM[1] = psiM[2];
            dPsiM[1] = dPsiM[2];
            psiM[2] = RK(xm, psiM[1], dPsiM[1], hX, E).first;
            dPsiM[2] = RK(xm, psiM[1], dPsiM[1], hX, E).second;
            xm = xRight + hX * k;
        }

        // Если значение энергии является особственным функции будут сшиваться
        // Находим определитель Вронского
        W[i] = dPsiN[2] * psiM[2] - psiN[2] * dPsiM[2];
        fout1 << E << " " << W[i] << endl;
    }
    fout1.close();

    // Отбор энергий
    double d = 0.0;
    int n = 0;
    for (int i = 1; i <= M; ++i) {
        if (W[i] * W[i - 1] <= 0) {
            d = d + hEnergy;
            Energy[n] = d;
            n++;
        } else d = d + hEnergy;
    }

    ///////////////////
    // Функция слева //
    ///////////////////
    n = 0;
    E = Energy[n];
    xRight = -xLeft;

    b = x0;
    a = xLeft;
    hX = (b - a) / N;

    // Волновая функция и ее производная на асимптотиках
    psiLeft[1] = exp(xLeft * sqrt(2.0 * E));
    dPsiLeft[1] = sqrt(2.0 * E) * exp(xLeft * sqrt(2.0 * E));
    xn = xLeft;

    for (int k = 0; k < N; ++k) {
        psiLeft[k + 1] = RK(xn, psiLeft[k], dPsiLeft[k], hX, E).first;
        dPsiLeft[k + 1] = RK(xn, psiLeft[k], dPsiLeft[k], hX, E).second;
        xn = xLeft + hX * k;
    }

    psiLeft[N + 1] = RK(xn, psiLeft[N], dPsiLeft[N], hX, E).first;

    S1 = pow(psiLeft[1], 2) - pow(psiLeft[N + 1], 2);
    for (int i = 2; i <= N; i += 2) {
        S1 = S1 + 4.0 * pow(psiLeft[i], 2) + 2.0 * pow(psiLeft[i + 1], 2);
    }

    S1 = S1 * hX / 3.0;
    //cout << S1 << endl;

    ofstream fout2;
    fout2.open(R"(..\\..\\NumericalMethods\\Schrodinger\\PsiLeft.csv)");
    for (int i = 1; i <= N + 1; ++i) {
        psiLeft[i] = pow(psiLeft[i], 2) / (S1 + S2);
        xn = xLeft + hX * (i - 1);
        fout2 << xn << " " << psiLeft[i] << endl;
    }
    fout2.close();

    ////////////////////
    // Функция справа //
    ////////////////////

    b = xRight;
    a = x0;
    h2 = (a - b) / N;

    // Волновая функция и ее производная на асимптотиках
    psiRight[1] = exp(-xRight * sqrt(2.0 * E));
    dPsiRight[1] = -sqrt(2.0 * E) * exp(-xRight * sqrt(2.0 * E));
    xm = xRight;

    for (int k = 1; k < N; ++k) {
        psiRight[k + 1] = RK(xm, psiRight[k], dPsiRight[k], h2, E).first;
        dPsiRight[k + 1] = RK(xm, psiRight[k], dPsiRight[k], h2, E).second;
        xm = xRight + h2 * k;
    }

    psiRight[N + 1] = RK(xm, psiRight[N], dPsiRight[N], h2, E).first;

    S2 = -pow(psiRight[1], 2) + pow(psiRight[N + 1], 2);
    for (int i = 2; i <= N; ++i) {
        S2 = S2 + 4.0 * pow(psiRight[i], 2) + 2.0 * pow(psiRight[i + 1], 2);
    }

    S2 = -1 * S2 * h2 / 3.0;
    //cout << S2 << endl;

    ofstream fout3;
    fout3.open(R"(..\\..\\NumericalMethods\\Schrodinger\\PsiRight.csv)");
    for (int i = 1; i <= N + 1; ++i) {
        psiRight[i] = pow(psiRight[i], 2) / (S1 + S2);
        xm = xRight + h2 * (i - 1);
        fout3 << xm << " " << psiRight[i] << endl;
    }
    fout3.close();
}