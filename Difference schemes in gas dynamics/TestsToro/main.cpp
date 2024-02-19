#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <fstream>

// для каждого теста нужно написать свою функцию чтения файлов

/* ПАРАМЕТРЫ */
struct PVR{
    double p;
    double v;
    double r;
};

/* ВХОДНЫЕ ДАННЫЕ */
struct Parameters{
    int cells_number;           /* число ячеек */
    double stop_time;           /* момент времени, для которого строится точное решение */
    std::vector<PVR> params;    /* вектор переменных различных областей */
    std::vector<double> l;      /* вектор разделителей областей */
    double g;                   /* показатель адибаты */
    double CFL;                 /* число Куранта */
    int S_type;                 /* 1 - G, 2 - GK, 3 - GKR */
    double a;                   /* начало расчетной области */
    double b;                   /* конец расчетной области */
    int fo;                     /* частота вывода файлов */
    int q;                      /* 1 - с искуственной вязкостью, 0 - без искусственной вязкости */
    int bound; /* 0 - жесткая стенка, 1 - нет жесткой стенки */
};

struct Test{
    std::string name;
    double rL;
    double vL;
    double pL;
    double rR;
    double vR;
    double pR;
    int bound;
    double time;
};

struct IndividualProblem{
    double time = 0.25;
    int bound = 0;
    double e = 1.0e6;

    double x_start = -100.0;
    double rL = 1.0;
    double vL = 0.0;
    double pL = 1.0e5;

    double x_L; // -le/2
    double x_mid = 0.0;
    double x_R; // +le/2
    double rExp = 1.0;
    double vExp;
    double pExp = 1.0e5;
   

    double rR = 1.0;
    double vR = 0.0;
    double pR = 1.0e5;
    double x_end = 100.0;
};

/* ВЫБОР ТЕСТА */
Test test_chooser(int test_number) {
    /* НАЧАЛЬНЫЕ УСЛОВИЯ */
    Test test;

    switch (test_number)
    {
    case 1:
        test.name = "Sod";
        test.rL = 1.0;
        test.vL = 0.0;
        test.pL = 1.0;
        test.rR = 0.125;
        test.vR = 0.0;
        test.pR = 0.1;
        test.bound = 1;
        test.time = 0.25;
        break;

    case 2:
        test.name = "123";
        test.rL = 1.0;
        test.vL = -2.0;
        test.pL = 0.4;
        test.rR = 1.0;
        test.vR = 2.0;
        test.pR = 0.4;
        test.bound = 1;
        test.time = 0.15;
        break;

    case 3:
        test.name = "WoodwardColellaLeft";
        test.rL = 1.0;
        test.vL = 0.0;
        test.pL = 1000.0;
        test.rR = 1.0;
        test.vR = 0.0;
        test.pR = 0.01;
        test.bound = 1;
        test.time = 0.012;
        break;

    case 4:
        test.name = "WoodwardColellaRight";
        test.rL = 1.0;
        test.vL = 0.0;
        test.pL = 0.01;
        test.rR = 1.0;
        test.vR = 0.0;
        test.pR = 100.0;
        test.bound = 1;
        test.time = 0.035;
        break;

    case 5:
        test.name = "Combined";
        test.rL = 5.99924;
        test.vL = 19.5975;
        test.pL = 460.894;
        test.rR = 5.99242;
        test.vR = -6.19633;
        test.pR = 46.0950;
        test.bound = 1;
        test.time = 0.035;
        break;

    }

    return test;
}

/* ЧТЕНИЕ ВХОДНЫХ ДАННЫХ */
/* исходный код */
//void read_parameters(struct Parameters* params) {
void read_parameters(struct Parameters* params, Test test) {
    params->g = 1.4;
    /* границы областей расчета */
    params->a = 0.0;
    params->b = 1.0;
    
    params->l.push_back(params->a); // Левый край
    /* НАЧАЛЬНЫЕ УСЛОВИЯ */
    PVR pvr;
    /* область слева */
    pvr.p = test.pL;
    pvr.v = test.vL;
    pvr.r = test.rL;
    params->params.push_back(pvr);
    
    params->l.push_back(0.5); // Граница областей
    
    /* облаcть два */
    pvr.p = test.pR;
    pvr.v = test.vR;
    pvr.r = test.rR;
    
    params->params.push_back(pvr);
    
    params->l.push_back(params->b); // Правый край

    /* ПАРАМЕТРЫ РАСЧЕТА*/
    params->cells_number = 500;
    params->stop_time = test.time;
    params->CFL = 0.2;
    params->S_type = 1;
    params->fo = 50; // частота вывода файлов, выводится каждый 50-й шаг
    params->q = 0; //искусственная вязкость
    params->bound = test.bound; // тип границы 0 - жесткая стенка, 1 - нет жесткой стенки
}


/* БЛОК РАБОТЫ С ПАМЯТЬЮ */
/* Выделение памяти под одномерный массив элементов */
template <typename T>
void Allocate(int size, T * &mass){
    mass = new T[size];
}
template <class First, class... Other>
void Allocate(int size, First &first, Other&... other){
    Allocate(size, first);
    Allocate(size, other...);
}
/* Очищение памяти одномерного массива элементов */
template<class T>
void Delete(T* &mass){
    delete [] mass;
}
template <class First, class... Other>
void Delete(First &first, Other&... other){
    Delete(first);
    Delete(other...);
}

/* ПЕРЕХОД ОТ КОНСЕРВАТИВНЫХ ПЕРЕМЕННЫХ И ОБРАТНО */
void Convert_cons_to_noncons(struct Parameters *params, double& p, double& v, double& r, double& m, double& imp, double& e){
    double g = params->g;
    p = (g - 1.0) * (e - 0.5 * pow(imp, 2.0)/m);
    v = imp/m;
    r = m;
}
void Convert_noncons_to_cons(struct Parameters *params, double& p, double& v, double& r, double& m, double& imp, double& e){
    double g = params->g;
    m = r;
    imp = r * v;
    e = 0.5 * r * pow(v, 2.0) + p / (g - 1.0);
}

/* РАСЧЕТ СКОРОСТИ ЗВУКА */
double Sound_velocity(struct Parameters *params, double p, double r){
    double g = params->g;
    return std::sqrt(g * p / r);
}

/* ИНИЦИАЛИЗАЦИЯ СЕТКИ */
void Build_grid(struct Parameters *params, double *xc, double *x){
    double h;   /* шаг сетки */
    double right_boundary_x = params->b;
    double left_boundary_x = params->a;
    int cells_num = params->cells_number;
    h = (right_boundary_x - left_boundary_x ) / cells_num;
    /* координаты узлов */
    for(int i = 0; i < cells_num + 1; ++i){
        x[i] = left_boundary_x + i * h;
    }
    /* координаты центров ячеек */
    for(int i = 0; i < cells_num; ++i){
        xc[i] = 0.5 * (x[i] + x[i + 1]);
    }
}

void Diff_flux_ncons(struct Parameters *params, double p, double v, double r, double& Fm, double& Fimp, double& Fe){
    double m, imp, e; /* консервативные переменные */
    Convert_noncons_to_cons(params, p, v, r, m, imp, e);
    Fm = r*v;
    Fimp = Fm*v + p;
    Fe = (p + e) * v;
}

/* Расчет функции F, определяющей скорость газа на контактном разрыве, и ее производной по давлению среды DF
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   1999. - P. 158. - Subroutine PREFUN.
   + поправка на двучленное уравнение состояния: */
void Calc_F_and_DF(struct Parameters *params, double curr_press, double p, double v, double r, double c, double *F, double *DF){
    double g = params->g;           /* показатель адиабаты */
    double p_ratio, fg, q;          /* вспомогательные переменные */
    p_ratio = curr_press / p;
    if (curr_press <= p) {
        /* волна разрежения */
        fg = 2.0 / ( g - 1.0 );
        *F = fg * c * ( pow( p_ratio, 1.0 / fg / g ) - 1.0 );
        *DF = ( 1.0 / r / c ) * pow( p_ratio, - 0.5 * ( g + 1.0 ) / g );
    }
    else {
        /* ударная волна */
        q = sqrt( 0.5 * ( g + 1.0 ) / g * p_ratio + 0.5 * ( g - 1.0 ) / g );
        *F = ( curr_press - p ) / c / r / q;
        *DF = 0.25 * ( ( g + 1.0 ) * p_ratio + 3 * g - 1.0 ) / g / r / c / pow( q, 3.0 );
    }
}

/* Определение начального приближения для расчета давления на контактном разрыве
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 157. - Subroutine GUESSP.
   Возвращает искомое начальное приближения */
double Pressure_initial_guess(struct Parameters *params, double pl, double vl, double rl,
                              double cl, double pr, double vr, double rr, double cr){
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
    p_lin = std::max( 0.0, 0.5 * ( pl + pr ) - 0.125 * ( vr - vl ) * ( rl + rr ) * ( cl + cr ) );
    p_min = std::min( pl, pr );
    p_max = std::max( pl, pr );
    p_ratio = p_max / p_min;

    if ( ( p_ratio <= p_ratio_max ) &&
       ( ( p_min < p_lin && p_lin < p_max ) || ( fabs( p_min - p_lin ) < eps || fabs( p_max - p_lin ) < eps ) ) ) {
        /* Начальное приближение из линеаризованной задачи */
        return p_lin;
    } else {
        if ( p_lin < p_min ) {
            /* Начальное приближение по двум волнам разрежения
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 302. - Formula (9.32) + поправка на двучленное уравнение состояния */
            g1 = 0.5 * ( g - 1.0 ) / g;
            return pow( ( ( cl + cr - 0.5 * ( g - 1.0 ) * ( vr - vl ) ) / ( cl / pow( pl, g1 ) + cr / pow( pr, g1 ) ) ), 1.0 / g1 );
        } else {
            /* Начальное приближение по двум ударным волнам
               Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
               1999. - P. 128. - Formula (4.48) + поправка на двучленное уравнение состояния */
            g1 = 2.0 / ( g + 1.0 );
            g2 = ( g - 1.0 ) / ( g + 1.0 );
            p1 = std::sqrt( g1 / rl / ( g2 * pl + p_lin ) );
            p2 = std::sqrt( g1 / rr / ( g2 * pr + p_lin ) );
            return ( p1 * pl + p2 * pr - ( vr - vl ) ) / ( p1 + p2 );
        }
    }
}

/* Итерационная процедура расчета давления и скорости на контактном разрыве
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 3rd Edition. - Springer,
   2009. - P. 155. - Subroutine STARPU.*/

void Contact_pressure_velocity(struct Parameters *params, double pl, double vl, double rl, double cl,
                               double pr, double vr, double rr, double cr, double& p_cont, double& v_cont){
    double p_old;         /* значение давления на предыдущей итерации */
    double fl, fr;        /* значения функций */
    double fld, frd;      /* значения производных */
    int iter_num = 0;     /* количество проведенных итераций */
    int iter_max = 300;   /* максимальное количество итераций */
    double criteria;      /* переменная для определения сходимости */
    double g = params->g; /* показатель адиабаты */
    double eps = 1.e-8;
    if (2.0 * ( cl + cr ) / ( g - 1.0 ) <= vr - vl){
        /* случай возникновения вакуума */
        printf( "\nContact_pressure_velocity -> vacuum is generated\n" );
    }
    /* расчет начального приближения для давления */
    p_old = Pressure_initial_guess(params, pl, vl, rl, cl, pr, vr, rr, cr);
    if (p_old < 0.0){
        printf( "\nContact_pressure_velocity -> initial pressure guess is negative " );
    }
    /* решение нелинейного уравнения для нахождения давления на контактном разрыве методом Ньютона-Рафсона */
    do {
        Calc_F_and_DF(params, p_old, pl, vl, rl, cl, &fl, &fld);
        Calc_F_and_DF(params, p_old, pr, vr, rr, cr, &fr, &frd);
        p_cont = p_old - ( fl + fr + vr - vl ) / ( fld + frd );
        criteria = 2.0 * std::fabs( ( p_cont - p_old ) / ( p_cont + p_old ) );
        iter_num++;
        if (iter_num > iter_max){
            printf( "\nContact_pressure_velocity -> number of iterations exceeds the maximum value.\n" );
        }
        if (p_cont < 0.0){
            printf( "\nContact_pressure_velocity -> pressure is negative.\n" );
        }
        p_old = p_cont;
    } while (criteria > eps);
    /* скорость контактного разрыва */
    v_cont = 0.5 * ( vl + vr + fr - fl );
}


/* Функция отбора решения
   Код: Toro E.F. Riemann Solvers and Numerical Methods for Fluid Dynamics. - 2nd Edition. - Springer,
   1999. - P. 158. - Subroutine SAMPLE.
*/
void Sample_solid_solution(struct Parameters *params, double pl, double vl, double rl, double cl, double pr, double vr, double rr, double cr,
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
    g1 = 0.5 * ( params->g - 1.0 ) / params->g;
    g2 = 0.5 * ( params->g + 1.0 ) / params->g;
    g3 = 2.0 * params->g / ( params->g - 1.0 );
    g4 = 2.0 / ( params->g - 1.0 );
    g5 = 2.0 / ( params->g + 1.0 );
    g6 = ( params->g - 1.0 ) / ( params->g + 1.0 );
    g7 = 0.5 * ( params->g - 1.0 );

    if ( s <= v_cont ) {
        /* рассматриваемая точка - слева от контактного разрыва */
        if ( p_cont <= pl ) {
            /* левая волна разрежения */
            shl = vl - cl;
            if ( s <= shl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                cml = cl * pow( p_cont / pl, g1 );
                stl = v_cont - cml;
                if ( s > stl ) {
                    /* параметры слева от контактного разрыва */
                    r = rl * pow( p_cont / pl, 1.0 / params->g );
                    v = v_cont;
                    p = p_cont;
                }
                else {
                    /* параметры внутри левой волны разрежения */
                    v = g5 * ( cl + g7 * vl + s );
                    c = g5 * ( cl + g7 * ( vl - s ) );
                    r = rl * pow( c / cl, g4 );
                    p = pl * pow( c / cl, g3 );
                }
            }
        }
        else {
            /* левая ударная волна */
            p_ratio = p_cont / pl;
            sl = vl - cl * std::sqrt( g2 * p_ratio + g1 );
            if ( s <= sl ) {
                /* параметры слева от разрыва */
                r = rl;
                v = vl;
                p = pl;
            }
            else {
                /* параметры за левой ударной волной */
                r = rl * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
    }
    else {
        /* рассматриваемая точка - справа от контактного разрыва */
        if ( p_cont > pr ) {
            /* правая ударная волна */
            p_ratio = p_cont / pr;
            sr = vr + cr * std::sqrt( g2 * p_ratio + g1 );
            if ( s >= sr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
                /* параметры за правой ударной волной */
                r = rr * ( p_ratio + g6 ) / ( p_ratio * g6 + 1.0 );
                v = v_cont;
                p = p_cont;
            }
        }
        else {
            /* правая волна разрежения */
            shr = vr + cr;
            if ( s >= shr ) {
                /* параметры справа от разрыва */
                r = rr;
                v = vr;
                p = pr;
            }
            else {
               cmr = cr * pow( p_cont / pr, g1 );
               str = v_cont + cmr;
               if ( s <= str ) {
                   /* параметры справа от контактного разрыва */
                   r = rr * pow( p_cont / pr, 1.0 / params->g );
                   v = v_cont;
                   p = p_cont;
               }
               else {
                    /* параметры внутри правой волны разрежения */
                    v = g5 * ( - cr + g7 * vr + s );
                    c = g5 * ( cr - g7 * ( vr - s ) );
                    r = rr * pow( c / cr, g4 );
                    p = pr * pow( c / cr, g3 );
               }
            }
        }
    }
    /* формирование выходного вектора с результатом */
    r_res = r;
    v_res = v;
    p_res = p;
    
}

void Riman_solver(struct Parameters *params, double ml, double impl, double el, double mr, double impr, double er, double &p_res, double &v_res, double &r_res ){
    double cr, cl;
    double p_cont, v_cont;
    double p, v, r;
    double pl, vl, rl, pr, vr, rr;
    double g = params->g;

    Convert_cons_to_noncons(params, pl, vl, rl, ml, impl, el);
    Convert_cons_to_noncons(params, pr, vr, rr, mr, impr, er);
    cl = Sound_velocity(params, pl, rl);
    cr = Sound_velocity(params, pr, rr);

    if (2.0 * ( cl + cr ) / ( g - 1.0 ) <= vr - vl){
        /* случай возникновения вакуума */
        printf( "\nContact_pressure_velocity -> vacuum is generated\n" );
    }
     /* итерационная процедура расчета давления и скорости газа на контактном разрыве*/
    Contact_pressure_velocity(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont);
    /* отбор решения */
    Sample_solid_solution(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont, 0.0, p, v, r);
    p_res = p;
    v_res = v;
    r_res = r;
}

void Godunov_flux(struct Parameters *params, double ml, double impl, double el, double mr, double impr, double er, double& Fm, double& Fimp, double& Fe){
    double p, v, r;
    /*решение задачи о распаде разрыва*/
    Riman_solver(params, ml, impl, el, mr, impr, er, p, v, r);
    /* расчет потока Годунова по вектору неконсервативных переменных*/
    Diff_flux_ncons(params, p, v, r, Fm, Fimp, Fe);
}

/*mib_mod с учетом знаков */
double min_mod(double a, double b){
    if( a * b < 0)
        return 0;
    if (fabs(a) < fabs(b))
        return a;
    else
        return -b;
}

void Kolgan_reconstraction(struct Parameters *params, double ml, double impl, double el, double mc, double impc, double ec,
                            double mr,  double impr, double er, double& p_kolgan_l, double& p_kolgan_r, double& v_kolgan_l,
                            double& v_kolgan_r, double& r_kolgan_l, double& r_kolgan_r){
    
    double cr, cl;
    double p_cont, v_cont;
    double p, v, r;
    double pl, vl, rl, pc, vc, rc, pr, vr, rr;
    double g = params->g;
    double a, b, min_mod_res;
    Convert_cons_to_noncons(params, pl, vl, rl, ml, impl, el);
    Convert_cons_to_noncons(params, pc, vc, rc, mc, impc, ec);
    Convert_cons_to_noncons(params, pr, vr, rr, mr, impr, er);
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


double calc_time_step(struct Parameters *params, double* x, double* m, double* imp, double* e, int time_step_number){
    double new_step = 1000000;
    double p, v, r, c;
    double c_step;
    double CFL = params->CFL;
    int N = params->cells_number;
    for(int i = 0; i < N; ++i){
        Convert_cons_to_noncons(params, p, v, r, m[i], imp[i], e[i]);
        c = Sound_velocity(params, p, r);
        c_step = CFL*(x[i + 1] - x[i]) / (std::fabs(v) + c);
        if (c_step < new_step){
            new_step = c_step;
        }
    }
    return new_step;
}

void Init_solution(struct Parameters *params, double* xc, double* p, double* v, double* r, double* m, double* imp, double* e){
    for(int i = 0; i < params->cells_number; ++i){
        for(int j = 0; j < params->l.size() - 1; ++j){
            if (params->l[j] < xc[i] && xc[i] <= params->l[j + 1]){
                p[i] = params->params[j].p;
                v[i] = params->params[j].v;
                r[i] = params->params[j].r;
            }
        }
        Convert_noncons_to_cons(params, p[i], v[i], r[i], m[i], imp[i], e[i]);
    }
}

void Analityc(struct Parameters *params, double time, int N, double l, double x, double ml, double impl, double el, double mr, double impr, double er, double& p, double& v, double& r){
    double cr, cl;
    double p_cont, v_cont, r_cont;
    double p1, v1, r1, p2, v2, r2, p3;
    double e, D, XRW1, XRW2, X_cont, XSW;
    double pl, vl, rl, pr, vr, rr;
    double g = params->g;

    Convert_cons_to_noncons(params, pl, vl, rl, ml, impl, el);
    Convert_cons_to_noncons(params, pr, vr, rr, mr, impr, er);
    cl = Sound_velocity(params, pl, rl);
    cr = Sound_velocity(params, pr, rr);

    if (2.0 * ( cl + cr ) / ( g - 1.0 ) <= vr - vl){
        /* случай возникновения вакуума */
        printf( "\nContact_pressure_velocity -> vacuum is generated\n" );
    }
     /* итерационная процедура расчета давления и скорости газа на контактном разрыве*/
    Contact_pressure_velocity(params, pl, vl, rl, cl, pr, vr, rr, cr, p_cont, v_cont);
    r_cont = rr*((g - 1)*pr + (g + 1)*p_cont)/((g + 1)*pr + (g - 1)*p_cont);
    r2 = rl*pow(p_cont/pl, 1/g);
    v1 = (2/(g + 1))*(cl - (l - x)/time);
    r1 = rl*pow((1 - (g - 1)/2*v1/cl),2/(g - 1));
    p1 = pl*pow((1 - (g - 1)/2*v1/cl),2*g/(g - 1));
    D = (p_cont - pr)/(rr*v_cont);
    XRW1 = l - cl*time;
    XRW2 = l - (cl - (g + 1)*v_cont/2)*time;
    X_cont = l + v_cont*time;
    XSW = l + D*time;
    if (x < XRW1){
        p = pl;
        v = vl;
        r = rl;
    } else if (x < XRW2){
        p = p1;
        v = v1;
        r = r1;
    } else if (x < X_cont){
        p = p_cont;
        v = v_cont;
        r = r2;
    } else if (x < XSW){
        p = p_cont;
        v = v_cont;
        r = r_cont;
    } else {
        p = pr;
        v = vr;
        r = rr;
    }
}


void Out_a(struct Parameters *params, double time, int N, double l, double* x, double* p, double* v, double* r){
    char Name_file[100];
    sprintf(Name_file, ".out/csv/Out-%f-.csv", time);
    std::ofstream SurfaceFile(Name_file, std::ios::app);
    SurfaceFile << "x;p;v;r;e;pa;va;ra;ea;\n";
    double ml, impl, el, mr, impr, er, pa, va, ra;
    for(int i = 0; i < N; ++i){
        Convert_noncons_to_cons(params, params->params[0].p, params->params[0].v, params->params[0].r, ml, impl, el);
        Convert_noncons_to_cons(params, params->params[1].p, params->params[1].v, params->params[1].r, mr, impr, er);
        Analityc(params, time, N, l, x[i], ml, impl, el, mr, impr, er, pa, va, ra);
        SurfaceFile << x[i] << ";" << p[i] << ";" << v[i] << ";" << r[i] << ";" << p[i]/((params->g - 1)*r[i]) << ";" << pa << ";" << va << ";" << ra << ";" << pa/((params->g - 1)*ra) << "\n";
    }
}

void Out(struct Parameters *params, double time, double* x, double* p, double* v, double* r){
    int N = params->cells_number;
    char Name_file[100];
    
    sprintf(Name_file, "./out/csv/Out-%f-.csv", time);
    std::ofstream SurfaceFile(Name_file, std::ios::app);
    SurfaceFile << "x;p;v;r;e;\n";
    for(int i = 0; i < N; ++i){
        SurfaceFile << x[i] << ";" << p[i] << ";" << v[i] << ";" << r[i] << ";" << p[i]/((params->g - 1)*r[i]) << "\n";
    }
}

/* ГРАНИЧНЫЕ УСЛОВИЯ */
void Boundary(double m, double imp, double e, double& mb, double& impb, double& eb, int b_type){
    /* стенка */
    if (b_type == 0){
        mb = m;
        impb = -imp;
        eb = e;
    } else  if(b_type == 1){ /* свободная */
        mb = m;
        impb = imp;
        eb = e;
    }
}

int main() {
    // Сохраним все значения le и best_x для дальнейшего использования
    std::vector<double> le_values;
    std::vector<double> best_x_values;

    struct Parameters params;   /* структура с параметрами вычислительного эксперимента  */
    double *xc;                 /* массив координат центров ячеек сетки */
    double *x;                  /* массив координат узлов сетки */
    double *p, *v, *r;          /* массивы неконс переменных */
    double *m, *imp, *e;        /* массивы конс переменных */
    double *m_next, *imp_next, *e_next;
    double mb, impb, eb;        /* граничные значения конс переменных */
    double FmL, FimpL, FeL, FmR, FimpR, FeR; /* потоки конс переменных */

    double pl, vl, rl, pr, vr, rr;
    double ml, impl, el, mr, impr, er;
    double dm, dimp, de, *m05, *imp05, *e05;

    double *p_kolgan_l, *p_kolgan_r, *v_kolgan_l, *v_kolgan_r, *r_kolgan_l, *r_kolgan_r;
    double *m_kolgan_l, *m_kolgan_r, *imp_kolgan_l, *imp_kolgan_r, *e_kolgan_l, *e_kolgan_r;

    int step_number = 0;
    double time = 0;
    double dt;
    
    /* 1 - "Sod", 2 - "123", 3 - "WoodwardColellaLeft",
    4 - "WoodwardColellaRight", 5 - "Combined" */
    int test_number = 5;
    
    /* считывание файла с параметрами задачи */
    Test test = test_chooser(test_number);
    read_parameters(&params, test);

    int N = params.cells_number;
    int S_type = params.S_type;
    int fo = params.fo;
    /* выделение памяти под массивы */
    Allocate(N, xc, p, v, r, m, imp, e, m_next, imp_next, e_next,
            p_kolgan_l, p_kolgan_r, v_kolgan_l, v_kolgan_r, r_kolgan_l, r_kolgan_r,
            m_kolgan_l, m_kolgan_r, imp_kolgan_l, imp_kolgan_r, e_kolgan_l, e_kolgan_r,
            m05, imp05, e05);
    Allocate(N + 1, x);
    /* определение координат центров ячеек сетки */
    Build_grid(&params, xc, x);
    /* начальные условия */
    Init_solution(&params, xc, p, v, r, m, imp, e);
    /* вывод в файл начального распределения*/
    // Out_a(&params, time, N, params.l[1], xc, p, v, r);
    Out(&params, time, xc, p, v, r);
    
    /* цикл по времени*/
    while (time < params.stop_time){
        if (params.q == 1){
            for(int i = 0; i < N; ++i){
                if (i != N - 1){
                    double c1 = 0.5;
                    double c2 = 0.0;
                    p[i] = p[i] + c1*r[i]*pow(v[i + 1] - v[i], 2.0) + c2*r[i]*Sound_velocity(&params, p[i], r[i])*fabs(v[i + 1] - v[i]);
                    /*p[i] = (p[i] + (params.g + 1)/4*r[i]*pow(v[i + 1] - v[i], 2.0)
                    + std::sqrt(pow((params.g + 1)/4*r[i]*pow(v[i + 1] - v[i], 2.0), 2.0)
                    + pow(r[i]*Sound_velocity(&params, p[i], r[i])*(v[i + 1] - v[i]), 2.0)));*/
                }else{
                    p[i] = p[i - 1];
                }
                Convert_noncons_to_cons(&params, p[i], v[i], r[i], m[i], imp[i], e[i]);
            }
        }
        dt = calc_time_step(&params, x, m, imp, e, step_number);
        if (S_type == 2 || S_type == 3){
            /* РЕКОНСТРУКЦИЯ ПО КОЛГАНУ */
            for(int i = 0; i < N; ++i){
                if(i == 0){
                    Boundary(m[0], imp[0], e[0], mb, impb, eb , test.bound);
                    Convert_cons_to_noncons(&params, p[0], v[0], r[0], m[0], imp[0], e[0]);
                    
                    Kolgan_reconstraction(&params, mb, impb, eb, m[i], imp[i], e[i], m[i + 1], imp[i + 1], e[i + 1],
                                            p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i], v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i]);
                    Convert_noncons_to_cons(&params, p_kolgan_l[i], v_kolgan_l[i], r_kolgan_l[i], m_kolgan_l[i], imp_kolgan_l[i], e_kolgan_l[i]);
                    Convert_noncons_to_cons(&params, p_kolgan_r[i], v_kolgan_r[i], r_kolgan_r[i], m_kolgan_r[i], imp_kolgan_r[i], e_kolgan_r[i]);
                }else if (i == N - 1){
                    Boundary(m[N - 1], imp[N - 1], e[N - 1], mb, impb, eb , test.bound);
                    Convert_cons_to_noncons(&params, p[N - 1], v[N - 1], r[N - 1], m[N - 1], imp[N - 1], e[N - 1]);
                    Kolgan_reconstraction(&params, m[i - 1], imp[i - 1], e[i - 1], m[i], imp[i], e[i], mb, impb, eb,
                                            p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i], v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i]);
                    Convert_noncons_to_cons(&params, p_kolgan_l[i], v_kolgan_l[i], r_kolgan_l[i], m_kolgan_l[i], imp_kolgan_l[i], e_kolgan_l[i]);
                    Convert_noncons_to_cons(&params, p_kolgan_r[i], v_kolgan_r[i], r_kolgan_r[i], m_kolgan_r[i], imp_kolgan_r[i], e_kolgan_r[i]);
                }else{
                    Kolgan_reconstraction(&params, m[i - 1], imp[i - 1], e[i - 1], m[i], imp[i], e[i], m[i + 1], imp[i + 1], e[i + 1],
                                            p_kolgan_l[i], p_kolgan_r[i], v_kolgan_l[i], v_kolgan_r[i], r_kolgan_l[i], r_kolgan_r[i]);
                    Convert_noncons_to_cons(&params, p_kolgan_l[i], v_kolgan_l[i], r_kolgan_l[i], m_kolgan_l[i], imp_kolgan_l[i], e_kolgan_l[i]);
                    Convert_noncons_to_cons(&params, p_kolgan_r[i], v_kolgan_r[i], r_kolgan_r[i], m_kolgan_r[i], imp_kolgan_r[i], e_kolgan_r[i]);
                }
            }
        }
        if (S_type == 3){
                /* РОДИОНОВ */
                for(int i = 0; i < N; ++i){
                    /* ШАГ ПРЕДИКТОР */
                    Diff_flux_ncons(&params, p_kolgan_l[i], v_kolgan_l[i], r_kolgan_l[i], FmL, FimpL, FeL);
                    Diff_flux_ncons(&params, p_kolgan_r[i], v_kolgan_r[i], r_kolgan_r[i], FmR, FimpR, FeR);

                    m05[i] = m[i] - dt*(FmR - FmL)/(xc[i + 1] - xc[i]);
                    imp05[i] = imp[i] - dt*(FimpR - FimpL)/(xc[i + 1] - xc[i]);
                    e05[i] = e[i] - dt*(FeR - FeL)/(xc[i + 1] - xc[i]);
                }
                /* ШАГ КОРРЕКТОР */
                for(int i = 0; i < N; ++i){
                    /* посчитанные приращения */
                    dm = m_kolgan_r[i] - m_kolgan_l[i];
                    dimp = imp_kolgan_r[i] - imp_kolgan_l[i];
                    de = e_kolgan_r[i] - e_kolgan_l[i];
                    /* расчет входных параметров для задачи Римана */
                    m_kolgan_r[i] = 0.5 * (m[i] + m05[i]) + dm;
                    imp_kolgan_r[i] = 0.5 * (imp[i] + imp05[i]) + dimp;
                    e_kolgan_r[i] = 0.5 * (e[i] + e05[i]) + de;

                    m_kolgan_l[i] = 0.5 * (m[i] + m05[i]) - dm;
                    imp_kolgan_l[i] = 0.5 * (imp[i] + imp05[i]) - dimp;
                    e_kolgan_l[i] = 0.5 * (e[i] + e05[i]) - de;
                }
            }
        for(int i = 0; i < N; ++i){
            /* расчет потока через левую грань ячейки */
            if (S_type == 1){
                /* ГОДУНОВ */
                if (i != 0){
                    ml = m[i - 1];
                    impl = imp[i - 1];
                    el = e[i - 1];
                    mr = m[i];
                    impr = imp[i];
                    er = e[i];
                } else {
                
            
                    Boundary(m[0], imp[0], e[0], mb, impb, eb, test.bound);
                    ml = mb;
                    impl = impb;
                    el = eb;
                    mr = m[i];
                    impr = imp[i];
                    er = e[i];
                    Convert_cons_to_noncons(&params, p[0], v[0], r[0], m[0], imp[0], e[0]);
                    
            
                    
                }
            } else if (S_type == 2 || S_type == 3){
                /* ГОДУНОВ - КОЛГАН И РОДИОНОВ*/
                if (i != 0){
                    ml = m_kolgan_r[i - 1];
                    impl = imp_kolgan_r[i - 1];
                    el = e_kolgan_r[i - 1];
                    mr = m_kolgan_l[i];
                    impr = imp_kolgan_l[i];
                    er = e_kolgan_l[i];
                } else if (i == 0){
                    Boundary(m_kolgan_r[0], imp_kolgan_r[0], e_kolgan_r[0], mb, impb, eb , test.bound);
                    ml = mb;
                    impl = impb;
                    el = eb;
                    mr = m_kolgan_l[i];
                    impr = imp_kolgan_l[i];
                    er = e_kolgan_l[i];
                }
            }
            Godunov_flux(&params, ml, impl, el, mr, impr, er, FmL, FimpL, FeL);

            /* расчет потока через правую грань ячейки */
            if (S_type == 1){
                /* ГОДУНОВ */
                if (i != N - 1){
                    ml = m[i];
                    impl = imp[i];
                    el = e[i];
                    mr = m[i + 1];
                    impr = imp[i + 1];
                    er = e[i + 1];
                } else {
                    Boundary(m[N - 1], imp[N - 1], e[N - 1], mb, impb, eb , test.bound);
                    ml = m[i];
                    impl = imp[i];
                    el = e[i];
                    mr = mb;
                    impr = impb;
                    er = eb;
                }
            } else if (S_type == 2 || S_type == 3){
                /* ГОДУНОВ - КОЛГАН  и РОДИОНОВ*/
                if (i != N - 1){
                    ml = m_kolgan_r[i];
                    impl = imp_kolgan_r[i];
                    el = e_kolgan_r[i];
                    mr = m_kolgan_l[i + 1];
                    impr = imp_kolgan_l[i + 1];
                    er = e_kolgan_l[i + 1];
                } else if (i == N - 1){
                    ml = m_kolgan_r[i];
                    impl = imp_kolgan_r[i];
                    el = e_kolgan_r[i];
                    Boundary(m_kolgan_l[N - 1], imp_kolgan_l[N - 1], e_kolgan_l[N - 1], mb, impb, eb , test.bound);
                    mr = mb;
                    impr = impb;
                    er = eb;
                }
            }
            Godunov_flux(&params, ml, impl, el, mr, impr, er, FmR, FimpR, FeR);
            m_next[i] = m[i] - dt*(FmR - FmL)/(xc[i + 1] - xc[i]);
            imp_next[i] = imp[i] - dt*(FimpR - FimpL)/(xc[i + 1] - xc[i]);
            e_next[i] = e[i] - dt*(FeR - FeL)/(xc[i + 1] - xc[i]);
        }
        for(int i = 0; i < N; i++){
            m[i] = m_next[i];
            imp[i] = imp_next[i];
            e[i] = e_next[i];
            Convert_cons_to_noncons(&params, p[i], v[i], r[i], m[i], imp[i], e[i]);
        }
        time += dt;
        step_number ++;

        // Находим размер области в которой давление превышает Pmax = 1.5e5
        double p_max = 1.5e5;
        double x_size = 0;
        int N = 500; // params->cells_number;

        /* запись в файл */
        if (step_number%fo == 0){
            // Out_a(&params, time, N, params.l[1], xc, p, v, r);
            Out(&params, time, xc, p, v, r);
        }
    }

    /* освобождение памяти */
    Delete(x, xc, p, v, r, m, imp, e, m_next, imp_next, e_next,
    p_kolgan_l, p_kolgan_r, v_kolgan_l, v_kolgan_r, r_kolgan_l, r_kolgan_r,
    m_kolgan_l, m_kolgan_r, imp_kolgan_l, imp_kolgan_r, e_kolgan_l, e_kolgan_r,
    m05, imp05, e05);

    return 0;

}
