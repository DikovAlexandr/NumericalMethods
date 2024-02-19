#define _CRT_SECURE_NO_WARNINGS

#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <fstream>

/* ПОСТАНОВКА ЗАДАЧИ */
void InitialDataVar7(double**& x, double**& p, double**& rho, double***& v, 
                                    double**& u, double**& e, int Nx, int Ny, double gamma,
                                     double R, double x_centr, double y_centr, double Lx, double Ly)
{
    double ls = 0.4*Lx;
    double P1 = 3 * 1E5;
    double rho1 = 3;
    double v1 = 0;
    double e1 = P1 / rho1 / (gamma - 1);

    double P0 = 1E5;
    double rho0 = 1;
    double e0 = P0 / rho0 / (gamma - 1);

    double rho2 = 1E5;
    double e2 = P0 / rho2 / (gamma - 1); 


    for (int i = 0; i <= Nx + 1; i++) {
        for (int j = 0; j <= Ny + 1; j++) {
            if (x[0][i] < ls){
                p[i][j] = P1;
                rho[i][j] = rho1;
                v[0][i][j] = v1;
                v[1][i][j] = v1;
                e[i][j] = e1 + (v[0][i][j] * v[0][i][j] + v[1][i][j] * v[1][i][j]) / 2.;
                u[i][j] = e1;
            }
            else{
                p[i][j] = P0;
                rho[i][j] = rho0;
                v[0][i][j] = v1;
                v[1][i][j] = v1;
                e[i][j] = e0 + (v[0][i][j] * v[0][i][j] + v[1][i][j] * v[1][i][j]) / 2.;
                u[i][j] = e0;
            }
            if (pow(x[0][i] - 0.65*Lx, 2) + pow(x[1][j] - 0.5*Ly, 2) <= pow(0.25*Ly, 2)){
                p[i][j] = P0;
                rho[i][j] = rho2;
                v[0][i][j] = v1;
                v[1][i][j] = v1;
                e[i][j] = e2 + (v[0][i][j] * v[0][i][j] + v[1][i][j] * v[1][i][j]) / 2.;
                u[i][j] = e2;
            }
        }
    }
}


/* ПАРАМЕТРЫ */
struct PVR {
    double p;
    double v;
    double r;
};

/* ВХОДНЫЕ ДАННЫЕ */
struct Parameters {
    int cells_number_x;           /* число ячеек x*/
    int cells_number_y;           /* число ячеек y*/
    double stop_time;           /* момент времени, для которого строится точное решение */
    std::vector<PVR> params;    /* вектор переменных различных областей */
    std::vector<double> l;      /* вектор разделителей областей */
    double g;                   /* показатель адибаты */
    double CFL;                 /* число Куранта */
    int S_type;                 /* 1 - G, 2 - GK */
    double ax;                   /* начало расчетной области x*/
    double bx;                   /* конец расчетной области x*/
    double ay;                   /* начало расчетной области y*/
    double by;                   /* конец расчетной области y*/
    int fo;                     /* частота вывода файлов */
    double type_b_down;         /* тип нижней границы */
    double type_b_up;           /* тип верхней границы */
    double type_b_left;         /* тип левой границы */
    double type_b_right;        /* тип правой границы */

    double gx;                  /* x компонента ускорения свободного падения */
    double gy;                   /* y компонента ускорения свободного падения */

};

/* ЧТЕНИЕ ВХОДНЫХ ДАННЫХ */
void read_parameters(struct Parameters* params) {
    params->g = 1.4;
    /* НАЧАЛЬНЫЕ УСЛОВИЯ */
    params->ax = 0.0;
    params->bx = 1;
    params->ay = 0.0;
    params->by = 1;
    /* ПАРАМЕТРЫ РАСЧЕТА*/
    params->cells_number_x = 100;
    params->cells_number_y = 100;
    params->stop_time = 2;
    params->CFL = 0.2;
    params->S_type = 1;
    params->fo = 10000;
    params->type_b_down = 0;
    params->type_b_up = 0;
    params->type_b_left = 0;
    params->type_b_right = 0;

    params->gx = 0;
    params->gy = -100;

}

/* БЛОК РАБОТЫ С ПАМЯТЬЮ */
/* Выделение памяти под одномерный массив элементов */
template <typename T>
void Allocate(int size, T*& mass) {
    mass = new T[size];
}
template <typename T>
void Allocate(int Nx, int Ny, T**& mass) {
    T *data = new T[Nx * Ny];
    mass = new T * [Nx];
    for (int i = 0; i < Nx; ++i) {
        mass[i] = data + i * Ny;
    }
}
template <class First, class... Other>
void Allocate(int Nx, int Ny, First& first, Other&... other) {
    Allocate(Nx, Ny, first);
    Allocate(Nx, Ny, other...);
}

/* Очищение памяти одномерного массива элементов */
template<class T>
void Delete(T*& mass) {
    delete[] mass;
}
template<class T>
void Delete(int Nx, T**& mass) {
    delete[] *mass;
    delete[] mass;
}
template <class First, class... Other>
void Delete(int Nx, First& first, Other&... other) {
    Delete(Nx, first);
    Delete(Nx, other...);
}
template <class First, class... Other>
void Delete(First& first, Other&... other) {
    Delete(first);
    Delete(other...);
}

/* ПЕРЕХОД ОТ КОНСЕРВАТИВНЫХ ПЕРЕМЕННЫХ И ОБРАТНО */
void Convert_cons_to_noncons(struct Parameters* params, double& p, double& vx, double& vy, double& r, double& m, double& impx, double& impy, double& e) {
    double g = params->g;
    p = (g - 1.0) * (e - 0.5 * (pow(impx, 2.0) + pow(impy, 2.0)) / m);
    vx = impx / m;
    vy = impy / m;
    r = m;
}
void Convert_noncons_to_cons(struct Parameters* params, double& p, double& vx, double& vy, double& r, double& m, double& impx, double& impy, double& e) {
    double g = params->g;
    m = r;
    impx = r * vx;
    impy = r * vy;
    e = 0.5 * r * (pow(vx, 2.0) + pow(vy, 2.0)) + p / (g - 1.0);
}

/* РАСЧЕТ СКОРОСТИ ЗВУКА */
double Sound_velocity(struct Parameters* params, double p, double r) {
    double g = params->g;
    return std::sqrt(g * p / r);
}

/* ИНИЦИАЛИЗАЦИЯ СЕТКИ */
void Build_grid(struct Parameters* params, double* xc, double* yc, double* x, double* y) {
    double hx, hy;   /* шаг сетки */
    double right_boundary_x = params->bx;
    double left_boundary_x = params->ax;
    double right_boundary_y = params->by;
    double left_boundary_y = params->ay;
    int Nx = params->cells_number_x;
    int Ny = params->cells_number_y;
    hx = (right_boundary_x - left_boundary_x) / Nx;
    hy = (right_boundary_y - left_boundary_y) / Ny;
    /* координаты узлов */
    for (int i = 0; i < Nx + 1; ++i) {
        x[i] = left_boundary_x + i * hx;
    }
    /* координаты центров ячеек */
    for (int i = 0; i < Nx; ++i) {
        xc[i] = 0.5 * (x[i] + x[i + 1]);
    }
    /* координаты узлов */
    for (int i = 0; i < Ny + 1; ++i) {
        y[i] = left_boundary_y + i * hy;
    }
    /* координаты центров ячеек */
    for (int i = 0; i < Ny; ++i) {
        yc[i] = 0.5 * (y[i] + y[i + 1]);
    }
}

/* ГРАНИЧНЫЕ УСЛОВИЯ */
void Boundary_x(double p, double vx, double vy, double r, double& pb, double& vxb, double& vyb, double& rb, int b_type) {
    /* стенка */
    if (b_type == 0) {
        pb = p;
        vxb = -vx;
        vyb = vy;
        rb = r;
    }
    else  if (b_type == 1) {/* свободная */
        pb = p;
        vxb = vx;
        vyb = vy;
        rb = r;
    }
}
void Boundary_y(double p, double vx, double vy, double r, double& pb, double& vxb, double& vyb, double& rb, int b_type) {
    /* стенка */
    if (b_type == 0) {
        pb = p;
        vxb = vx;
        vyb = -vy;
        rb = r;
    }
    else  if (b_type == 1) {/* свободная */
        pb = p;
        vxb = vx;
        vyb = vy;
        rb = r;
    }
}


/* Расчет функции F, определяющей скорость газа на контактном разрыве, и ее производной по давлению среды DF
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.
   + поправка на двучленное уравнение состояния: */
void Calc_F_and_DF(struct Parameters* params, double curr_press, double p, double v, double r, double c, double* F, double* DF) {
    double g = params->g;           /* показатель адиабаты */
    double p_ratio, fg, q;          /* вспомогательные переменные */
    p_ratio = curr_press / p;
    if (curr_press <= p) {
        /* волна разрежения */
        fg = 2.0 / (g - 1.0);
        *F = fg * c * (pow(p_ratio, 1.0 / fg / g) - 1.0);
        *DF = (1.0 / r / c) * pow(p_ratio, -0.5 * (g + 1.0) / g);
    }
    else {
        /* ударная волна */
        q = sqrt(0.5 * (g + 1.0) / g * p_ratio + 0.5 * (g - 1.0) / g);
        *F = (curr_press - p) / c / r / q;
        *DF = 0.25 * ((g + 1.0) * p_ratio + 3 * g - 1.0) / g / r / c / pow(q, 3.0);
    }
}

/* Определение начального приближения для расчета давления на контактном разрыве
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.
   Возвращает искомое начальное приближения */
double Pressure_initial_guess(struct Parameters* params, double pl, double vl, double rl,
    double cl, double pr, double vr, double rr, double cr) {
    double g = params->g;                           /* показатель адиабаты */

    /* начальное приближение, рассчитанное на освановании рассмотрения линеаризованной системы
       в примитивных переменных */
    double p_lin;
    double p_min, p_max;                /* минимальное и максимальное давления слева и справа от разрыва */
    double p_ratio;                     /* перепад по давлению слева и справа от разрыва */
    double p_ratio_max = 2.0;           /* максимальный перепад по давлению слева и справа от разрыва */
    double p1, p2, g1, g2;              /* вспомогательные переменные для промежуточных расчетов */
    double eps = 1.e-8;                 /* малый эпсилон для критерия сходимости итераций и сравнения вещественных чисел */

    /* Начальное приближение из линейной задачи
       Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
       1999. - P. 128. - Formula (4.47). */
    p_lin = std::max(0.0, 0.5 * (pl + pr) - 0.125 * (vr - vl) * (rl + rr) * (cl + cr));
    p_min = std::min(pl, pr);
    p_max = std::max(pl, pr);
    p_ratio = p_max / p_min;

    if ((p_ratio <= p_ratio_max) &&
        ((p_min < p_lin && p_lin < p_max) || (fabs(p_min - p_lin) < eps || fabs(p_max - p_lin) < eps))) {
        /* Начальное приближение из линеаризованной задачи */
        return p_lin;
    }
    else {
        if (p_lin < p_min) {
            /* Начальное приближение по двум волнам разрежения
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 302. - Formula (9.32) + поправка на двучленное уравнение состояния */
            g1 = 0.5 * (g - 1.0) / g;
            return pow(((cl + cr - 0.5 * (g - 1.0) * (vr - vl)) / (cl / pow(pl, g1) + cr / pow(pr, g1))), 1.0 / g1);
        }
        else {
            /* Начальное приближение по двум ударным волнам
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 128. - Formula (4.48) + поправка на двучленное уравнение состояния */
            g1 = 2.0 / (g + 1.0);
            g2 = (g - 1.0) / (g + 1.0);
            p1 = std::sqrt(g1 / rl / (g2 * pl + p_lin));
            p2 = std::sqrt(g1 / rr / (g2 * pr + p_lin));
            return (p1 * pl + p2 * pr - (vr - vl)) / (p1 + p2);
        }
    }

}

/* Итерационная процедура расчета давления и скорости на контактном разрыве
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   2009. - P. 155. - Subroutine STARPU.*/

void Contact_pressure_velocity(struct Parameters* params, double pl, double vl, double rl, double cl,
    double pr, double vr, double rr, double cr, double& p_cont, double& v_cont) {
    double p_old;         /* значение давления на предыдущей итерации */
    double fl, fr;        /* значения функций */
    double fld, frd;      /* значения производных */
    int iter_num = 0;     /* количество проведенных итераций */
    int iter_max = 300;   /* максимальное количество итераций */
    double criteria;      /* переменная для определения сходимости */
    double g = params->g; /* показатель адиабаты */
    double eps = 1.e-8;
    if (2.0 * (cl + cr) / (g - 1.0) <= vr - vl) {
        /* случай возникновения вакуума */
        printf("\nContact_pressure_velocity -> vacuum is generated\n");
    }
    /* расчет начального приближения для давления */
    p_old = Pressure_initial_guess(params, pl, vl, rl, cl, pr, vr, rr, cr);
    if (p_old < 0.0) {
        printf("\nContact_pressure_velocity -> initial pressure guess is negative ");
    }
    /* решение нелинейного уравнения для нахождения давления на контактном разрыве методом Ньютона-Рафсона */
    do {
        Calc_F_and_DF(params, p_old, pl, vl, rl, cl, &fl, &fld);
        Calc_F_and_DF(params, p_old, pr, vr, rr, cr, &fr, &frd);
        p_cont = p_old - (fl + fr + vr - vl) / (fld + frd);
        criteria = 2.0 * std::fabs((p_cont - p_old) / (p_cont + p_old));
        iter_num++;
        if (iter_num > iter_max) {
            printf("\nContact_pressure_velocity -> number of iterations exceeds the maximum value.\n");
        }
        if (p_cont < 0.0) {
            printf("\nContact_pressure_velocity -> pressure is negative.\n");
        }
        p_old = p_cont;
    } while (criteria > eps);
    /* скорость контактного разрыва */
    v_cont = 0.5 * (vl + vr + fr - fl);
}


/* Функция отбора решения
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
*/
void Sample_solid_solution(struct Parameters* params, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr,
    double p_cont, double v_cont, double s, double& p_res, double& v_res, double& r_res) {

    double g1, g2, g3, g4, g5, g6, g7;      /* вспомогательные переменные, производные от показателя адиабаты,
                                               в соответствии с Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer, 1999. - P. 153. */

                                               /* скорости левых волн */
    double shl, stl;        /* скорости "головы" и "хвоста" левой волны разрежения */
    double sl;              /* скорость левой ударной волны */

    /* скорости правых волн */
    double shr, str;        /* скорости "головы" и "хвоста" правой волны разрежения */
    double sr;              /* скорость правой ударной волны */

    double cml, cmr;        /* скорости звука слева и справа от контактного разрыва */
    double c;               /* локальная скорость звука внутри волны разрежения */
    double p_ratio;
    double r, v, p;         /* отобранные значения объемной доли, плотности, скорости и давления */

    /* производные от показателя адиабаты */
    g1 = 0.5 * (params->g - 1.0) / params->g;
    g2 = 0.5 * (params->g + 1.0) / params->g;
    g3 = 2.0 * params->g / (params->g - 1.0);
    g4 = 2.0 / (params->g - 1.0);
    g5 = 2.0 / (params->g + 1.0);
    g6 = (params->g - 1.0) / (params->g + 1.0);
    g7 = 0.5 * (params->g - 1.0);

    if (s <= v_cont) {
        /* рассматриваемая точка - слева от контактного разрыва */
        if (p_cont <= pl) {
            /* левая волна разрежения */
            shl = vl - cl;
            if (s <= shl) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow(p_cont / pl, g1);
                stl = v_cont - cml;
                if (s > stl) {
                    /* параметры слева от контактного разрыва */
                    r = rl * pow(p_cont / pl, 1.0 / params->g);
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* параметры внутри левой волны разрежения */
                    v = g5 * (cl + g7 * vl + s);
                    c = g5 * (cl + g7 * (vl - s));
                    r = rl * pow(c / cl, g4);
                    p = pl * pow(c / cl, g3);
                }
            }
        }
        else {
            /* левая ударная волна */
            p_ratio = p_cont / pl;
            sl = vl - cl * std::sqrt(g2 * p_ratio + g1);
            if (s <= sl) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* параметры за левой ударной волной */
                r = rl * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        /* рассматриваемая точка - справа от контактного разрыва */
        if (p_cont > pr) {
            /* правая ударная волна */
            p_ratio = p_cont / pr;
            sr = vr + cr * std::sqrt(g2 * p_ratio + g1);
            if (s >= sr) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* параметры за правой ударной волной */
                r = rr * (p_ratio + g6) / (p_ratio * g6 + 1.0);
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            /* правая волна разрежения */
            shr = vr + cr;
            if (s >= shr) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                cmr = cr * pow(p_cont / pr, g1);
                str = v_cont + cmr;
                if (s <= str) {
                    /* параметры справа от контактного разрыва */
                    r = rr * pow(p_cont / pr, 1.0 / params->g);
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* параметры внутри правой волны разрежения */
                    v = g5 * (-cr + g7 * vr + s);
                    c = g5 * (cr - g7 * (vr - s));
                    r = rr * pow(c / cr, g4);
                    p = pr * pow(c / cr, g3);
                }
            }
        }
    }
    /* формирование выходного вектора с результатом */
    r_res = r;
    v_res = v;
    p_res = p;

}

void Riman_solver(struct Parameters* params, double rl, double vl, double pl, double rr, double vr, double pr, double& p_res, double& v_res, double& r_res) {
    double cr, cl;
    double p_cont, v_cont;
    double p, v, r;
    double g = params->g;

    cl = Sound_velocity(params, pl, rl);
    cr = Sound_velocity(params, pr, rr);

    if (2.0 * (cl + cr) / (g - 1.0) <= vr - vl) {
        /* случай возникновения вакуума */
        printf("\nContact_pressure_velocity -> vacuum is generated\n");
    }
    /* итерационная процедура расчета давления и скорости газа на контактном разрыве*/
    Contact_pressure_velocity(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont);
    /* отбор решения */
    Sample_solid_solution(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont, 0.0, p, v, r);
    p_res = p;
    v_res = v;
    r_res = r;
}

void Diff_flux_ncons_x(struct Parameters* params, double p, double vx, double vy, double r, double& Fm, double& Fimp_x, double& Fimp_y, double& Fe) {
    double m, impx, impy, e; /* консервативные переменные */

    Convert_noncons_to_cons(params, p, vx, vy, r, m, impx, impy, e);
    Fm =  r * vx;
    Fimp_x = Fm * vx + p;
    Fimp_y = Fm * vy;
    Fe = (p + e) * vx;
}

void Diff_flux_ncons_y(struct Parameters* params, double p, double vx, double vy, double r, double& Fm, double& Fimp_x, double& Fimp_y, double& Fe) {
    double m, impx, impy, e; /* консервативные переменные */

    Convert_noncons_to_cons(params, p, vx, vy, r, m, impx, impy, e);
    Fm =  r * vy;
    Fimp_x = Fm * vx;
    Fimp_y = Fm * vy + p;
    Fe =  (p + e) * vy;
}

void Godunov_flux_x(struct Parameters* params, double ml, double impxl, double impyl, double el, double mr, double impxr, double impyr, double er, double& Fm, double& Fimp_x, double& Fimp_y, double& Fe) {
    double p, vx, vy, r;
    double pl, vxl, vyl, rl;
    double pr, vxr, vyr, rr;

    Convert_cons_to_noncons(params, pl, vxl, vyl, rl, ml, impxl, impyl, el);
    Convert_cons_to_noncons(params, pr, vxr, vyr, rr, mr, impxr, impyr, er);

    /*решение задачи о распаде разрыва*/
    Riman_solver(params, rl, vxl, pl, rr, vxr, pr, p, vx, r);

    if (vx >= 0)
        vy = impyl / ml;
    else
        vy = impyr / mr;

    /* расчет потока Годунова по вектору неконсервативных переменных*/
    Diff_flux_ncons_x(params, p, vx, vy, r, Fm, Fimp_x, Fimp_y, Fe);
}

void Godunov_flux_y(struct Parameters* params, double md, double impxd, double impyd, double ed, double mu, double impxu, double impyu, double eu, double& Fm, double& Fimp_x, double& Fimp_y, double& Fe) {
    double p, vx, vy, r;
    double pd, vxd, vyd, rd;
    double pu, vxu, vyu, ru;

    Convert_cons_to_noncons(params, pd, vxd, vyd, rd, md, impxd, impyd, ed);
    Convert_cons_to_noncons(params, pu, vxu, vyu, ru, mu, impxu, impyu, eu);

    /*решение задачи о распаде разрыва*/
    Riman_solver(params, rd, vyd, pd, ru, vyu, pu, p, vy, r);

    if (vy >= 0)
        vx = impxd / md;
    else
        vx = impxu / mu;

    /* расчет потока Годунова по вектору неконсервативных переменных*/
    Diff_flux_ncons_y(params, p, vx, vy, r, Fm, Fimp_x, Fimp_y, Fe);
}


/*mib_mod с учетом знаков */
double min_mod(double a, double b) {
    if (a * b < 0)
        return 0;
    if (fabs(a) < fabs(b))
        return a;
    else
        return -b;
}

// ##########################################################################################################################

void Kolgan_reconstraction(struct Parameters *params, 
                            double ml, double impxl, double impyl, double el, 
                            double mc, double impxc, double impyc, double ec, 
                            double mr, double impxr, double impyr, double er, 
                            double& p_kolgan_l, double& p_kolgan_r, double& v_kolgan_l, 
                            double& v_kolgan_r, double& r_kolgan_l, double& r_kolgan_r,
                            char direction){
    double pl, vxl, vyl, rl, pc, vxc, vyc, rc, pr, vxr, vyr, rr;
    // double g = params->g;
    double a, b, min_mod_res;
    // Convert_cons_to_noncons(params, p, vx, vy, r, m[i][j], impx[i][j], impy[i][j], e[i][j]);
    Convert_cons_to_noncons(params, pl, vxl, vyl, rl, ml, impxl, impyl, el);
    Convert_cons_to_noncons(params, pc, vxc, vyc, rc, mc, impxc, impyc, ec);
    Convert_cons_to_noncons(params, pr, vxr, vyr, rr, mr, impxr, impyr, er);
    
    double vl, impl, vc, impc, vr, impr;

    if (direction == 'x') {
        vl = vxl;
        impl = impxl;
        vc = vxc;
        impc = impxc;
        vr = vxr;
        impr = impxr;
    } 
    else {
        vl = vyl;
        impl = impyl;
        vc = vyc;
        impc = impyc;
        vr = vyr;
        impr = impyr;
    }
    /*плотность*/
    a = 0.5 * (rc - rl);
    b = 0.5 * (rr - rc);
    min_mod_res = min_mod((a + b)/2,  2*min_mod(a,b));
    r_kolgan_l = rc - min_mod_res;
    r_kolgan_r = rc + min_mod_res;
    /*давление*/
    a = 0.5 * (pc - pl);
    b = 0.5 * (pr - pc);
    min_mod_res = min_mod((a + b)/2,  2*min_mod(a,b));
    p_kolgan_l = pc - min_mod_res;
    p_kolgan_r = pc + min_mod_res;
    /*скорость*/
    a = 0.5 * (vc - vl);
    b = 0.5 * (vr - vc);
    min_mod_res = min_mod((a + b)/2,  2*min_mod(a,b));
    v_kolgan_l = vc - min_mod_res;
    v_kolgan_r = vc + min_mod_res;
}
// ##########################################################################################################################

double calc_time_step(struct Parameters* params, double* x, double* y, double** m, double** impx, double** impy, double** e, int time_step_number) {
    double new_step = 1000000;
    double p, vx, vy, r, c;
    double c_step;
    double CFL = params->CFL;
    for (int i = 0; i < params->cells_number_x; ++i) {
        for (int j = 0; j < params->cells_number_y; ++j) {
            Convert_cons_to_noncons(params, p, vx, vy, r, m[i][j], impx[i][j], impy[i][j], e[i][j]);
            c = Sound_velocity(params, p, r);
            c_step = std::min(CFL * (x[i + 1] - x[i]) / (std::fabs(vx) + c), CFL * (y[j + 1] - y[j]) / (std::fabs(vy) + c));
            if (c_step < new_step) {
                new_step = c_step;
            }
        }
    }
    return new_step;
}



void Init_solution_circle(struct Parameters* params, double* xc, double* yc, double** p, double** vx, double** vy, double** r, double** m, double** impx, double** impy, double** e) {
    double R = 0.1;
    double p1, vx1, vy1, r1, p2, vx2, vy2, r2;
   
    p1 = 1E5; // атмосферное давление
    vx1 = 0.;
    vy1 = 0;
    r1 = 0.175; // плотность воздуха

    p2 = 1E5; // атмосферное давление
    vx2 = 0;
    vy2 = 0;
    r2 = 0.1 * 0.175; // плотность воздуха


    for (int i = 0; i < params->cells_number_x; ++i) {
        for (int j = 0; j < params->cells_number_y; ++j) {
            if ((xc[i] - 0.5) * (xc[i] - 0.5) + (yc[j] - 0.2) * (yc[j] - 0.2) < R * R) {
            // if (yc[j] < R) {
            //if (xc[i] <  R) {
                p[i][j] = p1;
                vx[i][j] = vx1;
                vy[i][j] = vy1;
                r[i][j] = r1;
            }
            else {
                p[i][j] = p2;
                vx[i][j] = vx2;
                vy[i][j] = vy2;
                r[i][j] = r2;
            }
            Convert_noncons_to_cons(params, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
        }
    }
}


void Out(struct Parameters* params, double time, double* x, double* y, double** p, double** vx, double** vy, double** r) {
    int Nx = params->cells_number_x;
    int Ny = params->cells_number_y;
    char Name_file[100];
    sprintf(Name_file, "out/csv/Out_%f_.csv", time);
    std::ofstream SurfaceFile(Name_file, std::ios::app);
    SurfaceFile << "x;y;p;vx;vy;r;e;\n";
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            SurfaceFile << x[i] << ";" << y[j] << ";" << p[i][j] << ";" << vx[i][j] << ";" << vy[i][j] << ";" << r[i][j] << ";" << p[i][j] / ((params->g - 1) * r[i][j]) << "\n";
        }
    }
}


int main() {
    struct Parameters params;   /* структура с параметрами вычислительного эксперимента  */
    double* xc;                 /* массив x-координат центров ячеек сетки */
    double* x;                  /* массив x-координат узлов сетки */
    double* yc;                 /* массив y-координат центров ячеек сетки */
    double* y;                  /* массив y-координат узлов сетки */
    double** p, ** vx, ** vy, ** r;          /* массивы неконс переменных */
    double** m, ** impx, ** impy, ** e;        /* массивы конс переменных */
    double** m_next, ** impx_next, ** impy_next, ** e_next;
    double mb, impxb, impyb, eb;        /* граничные значения конс переменных */
    double FmL, FimpxL, FimpyL, FeL, FmR, FimpxR, FimpyR, FeR; /* потоки конс переменных */

    double pl, vl, rl, pr, vr, rr;
    double ml, impxl, impyl, el, mr, impxr, impyr, er;
    double dm, dimp, de;

    double *p_kolgan_l, *p_kolgan_r, *v_kolgan_l, *v_kolgan_r, *r_kolgan_l, *r_kolgan_r;  
    double *m_kolgan_l, *m_kolgan_r, *imp_kolgan_l, *imp_kolgan_r, *e_kolgan_l, *e_kolgan_r;

    int step_number = 0;
    double time = 0;
    double dt, dx, dy;
    /* считывание файла с параметрами задачи */
    read_parameters(&params);
    int Nx = params.cells_number_x;
    int Ny = params.cells_number_y;
    int S_type = params.S_type;
    int fo = params.fo;
    int type_b_down = params.type_b_down;
    int type_b_up = params.type_b_up;
    int type_b_left = params.type_b_left;
    int type_b_right = params.type_b_right;
    /* выделение памяти под массивы */
    Allocate(Nx, Ny, p, vx, vy, r, m, impx, impy, e, m_next, impx_next, impy_next, e_next); 
    Allocate(Nx + 1, xc);
    Allocate(Ny + 1, yc);
    Allocate(Nx + 1, x);
    Allocate(Ny + 1, y);

    /* определение координат центров ячеек сетки */
    Build_grid(&params, xc, yc, x, y);
    /* начальные условия */
    Init_solution_circle(&params, xc, yc, p, vx, vy, r, m, impx, impy, e);
    /* вывод в файл начального распределения*/
    Out(&params, time, xc, yc, p, vx, vy, r);
    /* цикл по времени*/
    while (time < params.stop_time) {
        dt = calc_time_step(&params, x, y, m, impx, impy, e, step_number);
        if (S_type == 2) {
            Allocate(Nx, p_kolgan_l);
            Allocate(Nx, p_kolgan_r);
            Allocate(Nx, v_kolgan_l);
            Allocate(Nx, v_kolgan_r);
            Allocate(Nx, r_kolgan_l);
            Allocate(Nx, r_kolgan_r);
            Allocate(Nx, m_kolgan_l);
            Allocate(Nx, m_kolgan_r);
            Allocate(Nx, imp_kolgan_l);
            Allocate(Nx, imp_kolgan_r);
            Allocate(Nx, e_kolgan_l);
            Allocate(Nx, e_kolgan_r);
            /* РЕКОНСТРУКЦИЯ ПО КОЛГАНУ вдоль x*/
            for(int i = 0; i < Nx; ++i){
                for(int j = 0; j < Ny; ++j){
                    if(i == 0){
                        Boundary_x(m[0][j], impx[0][j], impy[0][j], e[0][j], mb, impxb, impyb, eb, params.type_b_left);
                        Kolgan_reconstraction(&params, mb, impxb, impyb, eb,
                                                m[i][j], impx[i][j], impy[i][j], e[i][j],
                                                m[i + 1][j], impx[i + 1][j], impy[i + 1][j], e[i + 1][j],
                                                p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i],
                                                v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i],
                                                'x');
                        Convert_noncons_to_cons(&params, p_kolgan_l[i], v_kolgan_l[i], vy[i][j], r_kolgan_l[i], m_kolgan_l[i], imp_kolgan_l[i], impy[i][j], e_kolgan_l[i]);
                        Convert_noncons_to_cons(&params, p_kolgan_r[i], v_kolgan_r[i], vy[i][j], r_kolgan_r[i], m_kolgan_r[i], imp_kolgan_r[i], impy[i][j], e_kolgan_r[i]);
                    }else if (i == Nx - 1){
                        Boundary_x(m[Nx - 1][j], impx[Nx - 1][j], impy[Nx - 1][j], e[Nx - 1][j], mb, impxb, impyb, eb, params.type_b_right);
                        Kolgan_reconstraction(&params, m[i - 1][j], impx[i - 1][j], impy[i - 1][j], e[i - 1][j],
                                                m[i][j], impx[i][j], impy[i][j], e[i][j],
                                                mb, impxb, impyb, eb,
                                                p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i],
                                                v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i],
                                                'x');
                        Convert_noncons_to_cons(&params, p_kolgan_l[i], v_kolgan_l[i], vy[i][j], r_kolgan_l[i], m_kolgan_l[i], imp_kolgan_l[i], impy[i][j], e_kolgan_l[i]);
                        Convert_noncons_to_cons(&params, p_kolgan_r[i], v_kolgan_r[i], vy[i][j], r_kolgan_r[i], m_kolgan_r[i], imp_kolgan_r[i], impy[i][j], e_kolgan_r[i]);
                    }else{
                        Kolgan_reconstraction(&params, m[i - 1][j], impx[i - 1][j], impy[i - 1][j], e[i - 1][j],
                                                m[i][j], impx[i][j], impy[i][j], e[i][j],
                                                m[i + 1][j], impx[i + 1][j], impy[i + 1][j], e[i + 1][j],
                                                p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i], 
                                                v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i],
                                                'x');
                        Convert_noncons_to_cons(&params, p_kolgan_l[i], v_kolgan_l[i], vy[i][j], r_kolgan_l[i], m_kolgan_l[i], imp_kolgan_l[i], impy[i][j], e_kolgan_l[i]);
                        Convert_noncons_to_cons(&params, p_kolgan_r[i], v_kolgan_r[i], vy[i][j], r_kolgan_r[i], m_kolgan_r[i], imp_kolgan_r[i], impy[i][j], e_kolgan_r[i]);                        
                    }
                }
            }
            /* РЕКОНСТРУКЦИЯ ПО КОЛГАНУ вдоль y*/
            for(int j = 0; j < Ny; ++j){
                for(int i = 0; i < Nx; ++i){
                    if(j == 0){
                        Boundary_y(m[i][0], impx[i][0], impy[i][0], e[i][0], mb, impxb, impyb, eb, params.type_b_left);
                        Kolgan_reconstraction(&params, mb, impxb, impyb, eb,
                                                m[i][j], impx[i][j], impy[i][j], e[i][j],
                                                m[i][j + 1], impx[i][j + 1], impy[i][j + 1], e[i][j + 1],
                                                p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i],
                                                v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i],
                                                'y');
                        Convert_noncons_to_cons(&params, p_kolgan_l[i], vx[i][j], v_kolgan_l[i], r_kolgan_l[i], m_kolgan_l[i], impx[i][j], imp_kolgan_l[i], e_kolgan_l[i]);
                        Convert_noncons_to_cons(&params, p_kolgan_r[i], vx[i][j], v_kolgan_r[i], r_kolgan_r[i], m_kolgan_r[i], impx[i][j], imp_kolgan_r[i], e_kolgan_r[i]);
                    }else if (j == Ny - 1){
                        Boundary_y(m[i][Ny - 1], impx[i][Ny - 1], impy[i][Ny - 1], e[i][Ny - 1], mb, impxb, impyb, eb, params.type_b_right);
                        Kolgan_reconstraction(&params, m[i][j - 1], impx[i][j - 1], impy[i][j - 1], e[i][j - 1],
                                                m[i][j], impx[i][j], impy[i][j], e[i][j],
                                                mb, impxb, impyb, eb,
                                                p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i],
                                                v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i],
                                                'y');
                        Convert_noncons_to_cons(&params, p_kolgan_l[i], vx[i][j], v_kolgan_l[i], r_kolgan_l[i], m_kolgan_l[i], impx[i][j], imp_kolgan_l[i], e_kolgan_l[i]);
                        Convert_noncons_to_cons(&params, p_kolgan_r[i], vx[i][j], v_kolgan_r[i], r_kolgan_r[i], m_kolgan_r[i], impx[i][j], imp_kolgan_r[i], e_kolgan_r[i]);                        
                    }else{
                        Kolgan_reconstraction(&params, m[i][j - 1], impx[i][j - 1], impy[i][j - 1], e[i][j - 1],
                                                m[i][j], impx[i][j], impy[i][j], e[i][j],
                                                m[i][j + 1], impx[i][j + 1], impy[i][j + 1], e[i][j + 1],
                                                p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i], 
                                                v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i],
                                                'y');

                        Convert_noncons_to_cons(&params, p_kolgan_l[i], vx[i][j], v_kolgan_l[i], r_kolgan_l[i], m_kolgan_l[i], impx[i][j], imp_kolgan_l[i], e_kolgan_l[i]);
                        Convert_noncons_to_cons(&params, p_kolgan_r[i], vx[i][j], v_kolgan_r[i], r_kolgan_r[i], m_kolgan_r[i], impx[i][j], imp_kolgan_r[i], e_kolgan_r[i]);
                    }
                }
            }
            Delete(p_kolgan_l, p_kolgan_r, v_kolgan_l, v_kolgan_r, r_kolgan_l, r_kolgan_r, m_kolgan_l, m_kolgan_r, imp_kolgan_l, imp_kolgan_r, e_kolgan_l, e_kolgan_r);
        }
        for (int i = 0; i < Nx; ++i)
            for (int j = 0; j < Ny; ++j) {
                /* ПО ОСИ X */
                /* расчет потока через левую грань ячейки */
                if (S_type == 1) {
                    /* ГОДУНОВ */
                    if (i != 0) {
                        ml = m[i - 1][j];
                        impxl = impx[i - 1][j];
                        impyl = impy[i - 1][j];
                        el = e[i - 1][j];
                        mr = m[i][j];
                        impxr = impx[i][j];
                        impyr = impy[i][j];
                        er = e[i][j];
                    }
                    else {
                        Boundary_x(m[0][j], impx[0][j], impy[0][j], e[0][j], mb, impxb, impyb, eb, params.type_b_left);
                        ml = mb;
                        impxl = impxb;
                        impyl = impyb;
                        el = eb;
                        mr = m[i][j];
                        impxr = impx[i][j];
                        impyr = impy[i][j];
                        er = e[i][j];
                    }
                }
                else if (S_type == 2) {
                    /* ГОДУНОВ - КОЛГАН */
                }
                Godunov_flux_x(&params, ml, impxl, impyl, el, mr, impxr, impyr, er, FmL, FimpxL, FimpyL, FeL);

                /* расчет потока через правую грань ячейки */
                if (S_type == 1) {
                    /* ГОДУНОВ */
                    if (i != Nx - 1) {
                        ml = m[i][j];
                        impxl = impx[i][j];
                        impyl = impy[i][j];
                        el = e[i][j];
                        mr = m[i + 1][j];
                        impxr = impx[i + 1][j];
                        impyr = impy[i + 1][j];
                        er = e[i + 1][j];
                    }
                    else {
                        Boundary_x(m[Nx - 1][j], impx[Nx - 1][j], impy[Nx - 1][j], e[Nx - 1][j], mb, impxb, impyb, eb, params.type_b_right);
                        ml = m[i][j];
                        impxl = impx[i][j];
                        impyl = impy[i][j];
                        el = e[i][j];
                        mr = mb;
                        impxr = impxb;
                        impyr = impyb;
                        er = eb;
                    }
                }
                else if (S_type == 2){
                    /* ГОДУНОВ - КОЛГАН*/
                }
                Godunov_flux_x(&params, ml, impxl, impyl, el, mr, impxr, impyr, er, FmR, FimpxR, FimpyR, FeR);

                if (i != Nx - 1)
                    dx = (xc[i + 1] - xc[i]);
                else
                    dx = (xc[i] - xc[i - 1]);
                m_next[i][j] = m[i][j] - dt * (FmR - FmL) / dx;
                impx_next[i][j] = impx[i][j] - dt * (FimpxR - FimpxL) / dx;
                impy_next[i][j] = impy[i][j] - dt * (FimpyR - FimpyL) / dx;
                e_next[i][j] = e[i][j] - dt * (FeR - FeL) / dx;

                /* ПО ОСИ Y */
                /* расчет потока через нижнюю грань ячейки */
                if (S_type == 1) {
                    /* ГОДУНОВ */
                    if (j != 0) {
                        ml = m[i][j - 1];
                        impxl = impx[i][j - 1];
                        impyl = impy[i][j - 1];
                        el = e[i][j - 1];
                        mr = m[i][j];
                        impxr = impx[i][j];
                        impyr = impy[i][j];
                        er = e[i][j];
                    }
                    else {
                            Boundary_y(m[i][0], impx[i][0], impy[i][0], e[i][0], mb, impxb, impyb, eb, params.type_b_down);
                            ml = mb;
                            impxl = impxb;
                            impyl = impyb;
                            el = eb;
                            mr = m[i][j];
                            impxr = impx[i][j];
                            impyr = impy[i][j];
                            er = e[i][j];
                    }
                }
                else if (S_type == 2){
                    /* ГОДУНОВ - КОЛГАН */
                }
                Godunov_flux_y(&params, ml, impxl, impyl, el, mr, impxr, impyr, er, FmL, FimpxL, FimpyL, FeL);

                /* расчет потока через верхнюю грань ячейки */
                if (S_type == 1) {
                    /* ГОДУНОВ */
                    if (j != Ny - 1) {
                        ml = m[i][j];
                        impxl = impx[i][j];
                        impyl = impy[i][j];
                        el = e[i][j];
                        mr = m[i][j + 1];
                        impxr = impx[i][j + 1];
                        impyr = impy[i][j + 1];
                        er = e[i][j + 1];
                    }
                    else {
                        Boundary_y(m[i][Ny - 1], impx[i][Ny - 1], impy[i][Ny - 1], e[i][Ny - 1], mb, impxb, impyb, eb, params.type_b_up);
                        ml = m[i][j];
                        impxl = impx[i][j];
                        impyl = impy[i][j];
                        el = e[i][j];
                        mr = mb;
                        impxr = impxb;
                        impyr = impyb;
                        er = eb;
                    }
                }
                else if (S_type == 2) {
                    /* ГОДУНОВ - КОЛГАН */
                }
                Godunov_flux_y(&params, ml, impxl, impyl, el, mr, impxr, impyr, er, FmR, FimpxR, FimpyR, FeR);

                if (j != Ny - 1)
                    dy = (yc[j + 1] - yc[j]);
                else
                    dy = (yc[j] - yc[j - 1]);
                m_next[i][j] = m_next[i][j] - dt * (FmR - FmL) / dy;
                impx_next[i][j] = impx_next[i][j] - dt * (FimpxR - FimpxL) / dy;
                impy_next[i][j] = impy_next[i][j] - dt * (FimpyR - FimpyL) / dy;
                e_next[i][j] = e_next[i][j] - dt * (FeR - FeL) / dy;
            }


        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {               
                impx[i][j] = impx_next[i][j] + m[i][j] * dt * params.gx;
                impy[i][j] = impy_next[i][j] + m[i][j] * dt * params.gy;
                m[i][j] = m_next[i][j];
                e[i][j] = e_next[i][j];
                Convert_cons_to_noncons(&params, p[i][j], vx[i][j], vy[i][j], r[i][j], m[i][j], impx[i][j], impy[i][j], e[i][j]);
            }
        }
        time += dt;
        step_number++;
        /* запись в файл */
        if (step_number % fo == 0) {
        // if (time >= 0.119363) {
            Out(&params, time, xc, yc, p, vx, vy, r);
         //   break;
        }
    }
    /* освобождение памяти */
    Delete(x, xc, p, vx, vy, r, m, impx, impy, e, m_next, impx_next, impy_next, e_next);
    return 0;
}