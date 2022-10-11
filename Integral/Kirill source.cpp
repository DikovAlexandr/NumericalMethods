#include <iostream>
#include <cmath>

using namespace std;

double int_func(double);
double* grid_simpson(double, double, double);
double  integration_simpson(double* arr, double ,double ,double );
double* grid_trap(double, double, double);
double  integration_trap(double* arr, double, double, double);
double* grid_trivos(double, double, double);
double  integration_trivos(double* arr, double, double, double);



int main() {
	double a = 0, b = 3,n = 5000;
	double * x;
	x = grid_simpson(a, b, n);
	cout.precision(10);
	cout <<"simson_int"<< integration_simpson(x, n, a, b) << endl;
	x = grid_trap(a, b, n);
	cout << "trap_int" << integration_trap(x, n, a, b) << endl;
	x = grid_trivos(a, b, n);
	cout << "trivos_int" << integration_trivos(x, n, a, b) << endl;

	return 0;
}


double int_func(double x) {
	return x*x*x;
}


double *grid_simpson(double a, double b, double n) {
	double h;
	h = (b - a) / (2 * n); // n chetnoe
	double x_step[100000];
	x_step[0] = a;
	for (int i = 1; i < 2 * n + 1; i++) {
		x_step[i] = x_step[i - 1] + h;
	}
	return x_step;
}


double integration_simpson(double* arr, double n ,double a,double b) {
	double h,I_s;
	h = (b - a) / (2 * n);
	I_s = (h / 3) * (int_func(a) - int_func(b));
	for (int i = 1; i < 2 * n + 1; i++) {
		if (i % 2 != 0) {
			I_s = I_s + (h / 3) * 4 * int_func(arr[i]);
		}
		else {
			I_s = I_s + (h / 3) * 2 * int_func(arr[i]);
		}

	}

	return I_s;
}

double* grid_trap(double a, double b, double n) {
	double h;
	h = (b - a) / (n); // n chetnoe
	double x_step[100000];
	x_step[0] = a;
	for (int i = 1; i < n + 1; i++) {
		x_step[i] = x_step[i - 1] + h;
	}
	return x_step;
}


double integration_trap(double* arr, double n, double a, double b) {
	double h, I_s;
	h = (b - a) / ( n);
	I_s = (h / 2) * (int_func(a) - int_func(b));
	for (int i = 1; i <  n + 1; i++) {
		I_s = I_s + (h) * int_func(arr[i]);
	}

	return I_s;
}


double* grid_trivos(double a, double b, double n) {
	double h;
	h = (b - a) / (3 * n); // n chetnoe
	double x_step[100000];
	x_step[0] = a;
	for (int i = 1; i < 3 * n + 1; i++) {
		x_step[i] = x_step[i - 1] + h;
	}
	return x_step;
}


double integration_trivos(double* arr, double n, double a, double b) {
	double h, I_s;
	h = (b - a) / (3 * n);
	I_s = ((3 * h) / 8) * (int_func(a) - int_func(b));
	for (int i = 1; i < 3 * n + 1; i++) {
		if (i % 3 != 0) {
			I_s = I_s + ((3 * h) / 8) * 3 * int_func(arr[i]);
		}
		else {
			I_s = I_s + ((3 * h) / 8) * 2 * int_func(arr[i]);
		}
	}

	return I_s;
 }