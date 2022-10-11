#include <iostream>
#include <cmath>
using namespace std;

// Интегрируемая функция
double f(double x, double y)
{
    return pow(x,2)*pow(y,2);
}

// Область интегрирования
double limit(double x)
{
    return sqrt(1-pow(x,2))/4;
}

// Двойное интегрирование
double doubleIntegral(double h, double lx, double ux)
{
    int nx, ny;
    // z хранит значения функции
    // ax[] хранит интеграл wrt y
    // для всех x точки принимаются
    double z[50][50], ax[50], answer;
  
    // Считаем количество точек
    // в х и у интеграл
    nx = int((ux - lx) / h + 1);
    ny = int((limit(0.0) - (-limit(0.0))) / h + 1);
  
    // Расчитываем значения функции в узловых точках
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            z[i][j] = f(lx + i * h, limit(-1.0) + j * h);
        }
    }
  
    // Расчитываем значение интеграла методом симпсона
    for (int i = 0; i < nx; ++i) {
        ax[i] = 0;
        for (int j = 0; j < ny; ++j) {
            if (j == 0 || j == ny - 1)
                ax[i] += z[i][j];
            else if (j % 2 == 0)
                ax[i] += 2 * z[i][j];
            else
                ax[i] += 4 * z[i][j];
        }
        ax[i] *= (h / 3);
    }
  
    answer = 0;
  
    // Вычисление итогового значения интеграла
    // используем интеграл полученный в шаге выше
    for (int i = 0; i < nx; ++i) {
        if (i == 0 || i == nx - 1)
            answer += ax[i];
        else if (i % 2 == 0)
            answer += 2 * ax[i];
        else
            answer += 4 * ax[i];
    }
    answer *= (h / 3);
  
    return answer;
}

int main()
{
    // lx и ux верхний и нижний пределы интегрирования по x
    // ly и uy верхний и нижний пределы интегрирования по y
    // h шаг интегрирования по x
    // k шаг интегрирования по y
    double h, k, lx, ly, ux, uy;



    lx = -1, ux = 1, h = 0.1;

    cout << doubleIntegral(h, lx, ux) << endl;
    return 0;
}