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
    return psi1 * 2.0 * (eigenEnergy - 1.0 / cosh(x));
}

/*
double g(double x, double psi1, double psi2, double eigenEnergy) {
    return psi1 * 2.0 * (eigenEnergy - 1.0 / pow(cosh(x), 2));
}
*/

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
    const int N = 1000, M = 1000, M1 = 10;
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
    double psiLeft[N + 2][2], dPsiLeft[N + 1][2], psiRight[N + 2][2], dPsiRight[N + 1][2];
    double Energy[N+2];
    double S1[N+2], S2[N+2];

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
    double energyI = 0.0;
    int n = 0;
    for (int i = 1; i <= M; ++i) {
        if (W[i] * W[i - 1] <= 0) {

            energyI = energyI + hEnergy;
            Energy[n] = energyI;
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
        hX = (b - a) / N;

        // Волновая функция и ее производная на асимптотиках
        psiLeft[1][i] = exp(xLeft * sqrt(2.0 * Energy[i]));
        dPsiLeft[1][i] = sqrt(2.0 * Energy[n]) * exp(xLeft * sqrt(2.0 * Energy[i]));
        xn = xLeft;

        for (int k = 1; k < N; ++k) {
            psiLeft[k + 1][i] = RK(xn, psiLeft[k][i], dPsiLeft[k][i], hX, Energy[i]).first;
            dPsiLeft[k + 1][i] = RK(xn, psiLeft[k][i], dPsiLeft[k][i], hX, Energy[i]).second;
            xn = xLeft + hX * k;
        }

        psiLeft[N + 1][i] = RK(xn, psiLeft[N][i], dPsiLeft[N][i], hX, Energy[i]).first;

        // Нормировка методом Симпсона
        S1[i] = pow(psiLeft[1][i], 2) - pow(psiLeft[N + 1][i], 2);
        for (int j = 2; j <= N; j += 2) {
            S1[i] = S1[i] + 4.0 * pow(psiLeft[j][i], 2) + 2.0 * pow(psiLeft[j + 1][i], 2);
        }
        S1[i] = S1[i] * hX / 3.0;

        ////////////////////
        // Функция справа //
        ////////////////////

        b = xRight;
        a = x0;
        h2 = (a - b) / N;

        // Волновая функция и ее производная на асимптотиках
        psiRight[1][i] = exp(-xRight * sqrt(2.0 * Energy[i]));
        dPsiRight[1][i] = -sqrt(2.0 * Energy[n]) * exp(-xRight * sqrt(2.0 * Energy[i]));
        xm = xRight;

        for (int k = 1; k < N; ++k) {
            psiRight[k + 1][i] = RK(xm, psiRight[k][i], dPsiRight[k][i], h2, Energy[i]).first;
            dPsiRight[k + 1][i] = RK(xm, psiRight[k][i], dPsiRight[k][i], h2, Energy[i]).second;
            xm = xRight + h2 * k;
        }

        psiRight[N + 1][i] = RK(xm, psiRight[N][i], dPsiRight[N][i], h2, Energy[i]).first;

        // Нормировка методом Симпсона
        S2[i] = -pow(psiRight[1][i], 2) + pow(psiRight[N + 1][i], 2);
        for (int j = 2; j <= N; ++j) {
            S2[i] = S2[i] + 4.0 * pow(psiRight[j][i], 2) + 2.0 * pow(psiRight[j + 1][i], 2);
        }
        S2[i] = -1 * S2[i] * h2 / 3.0;

    }

    // Значения норированной волновой функции слева
    ofstream fout2;
    fout2.open(R"(..\\..\\SchrodingerOneDimension\\PsiLeft.csv)");
    for (int i = 1; i <= N + 1; ++i) {
        xn = xLeft + hX * (i - 1);
        fout2 << xn;
        for (int j = 0; j < n; ++j) {
            //psiLeft[i][j] = pow(psiLeft[i][j], 2) / (S1[j] + S2[j]);
            psiLeft[i][j] = psiLeft[i][j] / (S1[j] + S2[j]);
            fout2 << " " << psiLeft[i][j];
        }
        fout2 << endl;
    }
    fout2.close();

    // Значения норированной волновой функции справа
    ofstream fout3;
    fout3.open(R"(..\\..\\SchrodingerOneDimension\\PsiRight.csv)");
    for (int i = 1; i <= N + 1; ++i) {
        xm = xRight + h2 * (i - 1);
        fout3 << xm;
        for (int j = 0; j < n; ++j) {
            // psiRight[i][j] = pow(psiRight[i][j], 2) / (S1[j] + S2[j]);
            psiRight[i][j] = psiRight[i][j] / (S1[j] + S2[j]);
            fout3 << " " << psiRight[i][j];
        }
        fout3 << endl;
    }
    fout3.close();
}