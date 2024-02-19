#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace std;

// Тип float используется, поскольку имеет 24 разрядную мантиссу

// Просто суммирование слогаемых для демонстрации накопления ошибки
float dumbSum(float value, int k) {
    float sum = 0;
    float save = 0;
    for (int i = k; i > 0; i--) {
        save = sum;
        sum += value;
        // На определенном шаге суммирование перестает происходить, можем вывести номер этого шага
        // if (save == sum) {
        //     cout << "Summation was stopped at step " << k - i + 1 << endl;
        //    break;
        //}
    }
    return sum;
}

// Реализация метода Кэхэна для нахождения значения экспоненты
float kahan(float value, int k) {
    float current;
    float element;
    float error = 0.0;
    float save = -1.0;
    float sum = 0.0;
    while (k > 0) {
        save = sum;
        element = value - error;
        current = sum + element;
        error = current - sum - element;
        sum = current;
        // Если мы хотим узнать за сколько шагов набирается эта сумма, используем следующее условие выхода
        // if (sum == dumbSum(value, k)) {
        //    cout <<  "Kahan was stopped at step " << int(1.0e9 - k + 1) << endl;
        //    break;
        //}
        k--;
    }
    return sum;
}

int main() {
    float value = 1.0e-9;
    int k = 1.0e9;

    // Настройки точности вывода
    const auto digits = 8;
    cout << setw(digits);
    cout << fixed << setprecision(digits);

    // Нахождение значения прямым суммированием
    float dumbSumValue = dumbSum(value, k);
    // Нахождение значения Кэхоном
    float summExpValue = kahan(value, k);
    float analyticValue = 1;
    cout << "Summation result = " << dumbSumValue << endl;
    cout << "Kahan result = " << summExpValue << endl;
    cout << "Analytics result = " << analyticValue << endl;
    cout << "Delta summation = " << abs(dumbSumValue - analyticValue) << endl;
    cout << "Delta Kahan = " << abs(summExpValue - analyticValue) << endl;
    cout << "N = " << k << endl;
    return 0;
}
