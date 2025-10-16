#include "methods.h"

#include <cmath>       // fabs
#include <fstream>     // ofstream (CSV)
#include <iomanip>     // setprecision, fixed
#include <stdexcept>   // exceptions

/***********************
 * Исходная функция и производная
 ***********************/

// f(x) = x^3 - 2
double f(double x) {
    return x*x*x - 2.0;
}

// f'(x) = 3x^2
double df(double x) {
    return 3.0 * x * x;
}

/***********************
 * 1) Метод половинного деления (бисекции)
 ***********************/
MethodResultBisection solve_bisection(double a, double b, double eps, int Nmax) {
    if (a >= b) throw std::invalid_argument("bisection: a >= b");
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0.0) {
        // Если знаки одинаковые — не гарантируется наличие корня на отрезке
        throw std::invalid_argument("bisection: f(a) and f(b) must have opposite signs");
    }

    MethodResultBisection res{};
    res.eps = eps;

    for (int n = 0; n < Nmax; ++n) {
        double c  = 0.5 * (a + b);
        double fc = f(c);

        // Лог итерации: интервал, середина, значение функции
        res.rows.push_back({n, a, b, c, fc});

        // Критерии останова:
        // 1) |f(c)| < eps  ИЛИ
        // 2) половина длины интервала < eps
        if (std::fabs(fc) < eps || 0.5 * (b - a) < eps) {
            res.root       = c;
            res.f_at_root  = fc;
            res.iterations = n + 1;
            return res;
        }

        // Выбор подынтервала в зависимости от знака
        if (fa * fc < 0.0) {
            b  = c;
            fb = fc;
        } else {
            a  = c;
            fa = fc;
        }
    }

    // Если вышли по Nmax — возвращаем текущее лучшее
    double c = 0.5 * (a + b);
    res.root       = c;
    res.f_at_root  = f(c);
    res.iterations = Nmax;
    return res;
}

/***********************
 * 2) Метод простых итераций (фиксированной точки)
 ***********************/
MethodResultIter solve_fixed_point(double x0, double eps, int Nmax) {
    // Контрактная функция: phi(x) = 0.5*(x + 2/x^2),
    // в корне |phi'(xi)| = 0.5 < 1 => локальная сходимость.
    auto phi = [](double x) {
        return 0.5 * (x + 2.0 / (x * x));
    };

    MethodResultIter res{};
    res.eps = eps;

    double x = x0;
    for (int n = 0; n < Nmax; ++n) {
        // Страховка от деления на 0
        if (std::fabs(x) < 1e-14) x = 1e-6;

        double xn1   = phi(x);
        double delta = std::fabs(xn1 - x);
        double r     = std::fabs(f(xn1));

        res.rows.push_back({n, xn1, f(xn1), delta, r});

        if (delta < eps || r < eps) {
            res.root       = xn1;
            res.f_at_root  = f(xn1);
            res.iterations = n + 1;
            return res;
        }
        x = xn1;
    }

    res.root       = x;
    res.f_at_root  = f(x);
    res.iterations = Nmax;
    return res;
}

/***********************
 * 3) Метод Ньютона
 ***********************/
MethodResultIter solve_newton(double x0, double eps, int Nmax) {
    MethodResultIter res{};
    res.eps = eps;

    double x = x0;
    for (int n = 0; n < Nmax; ++n) {
        double y  = f(x);
        double dy = df(x);

        // Если производная слишком мала — остановим итерации
        if (std::fabs(dy) < 1e-14) {
            break;
        }

        double xn1   = x - y / dy;
        double delta = std::fabs(xn1 - x);
        double r     = std::fabs(f(xn1));

        res.rows.push_back({n, xn1, f(xn1), delta, r});

        if (delta < eps || r < eps) {
            res.root       = xn1;
            res.f_at_root  = f(xn1);
            res.iterations = n + 1;
            return res;
        }
        x = xn1;
    }

    res.root       = x;
    res.f_at_root  = f(x);
    res.iterations = (int)res.rows.size();
    return res;
}

/***********************
 * Экспорт CSV (итерационные таблицы)
 ***********************/

// маленькая проверка открытия
static void ensure_open(std::ofstream& out) {
    if (!out.good()) {
        throw std::runtime_error("Cannot open CSV file for writing");
    }
}

// CSV для бисекции
void save_bisection_csv(const std::string& path, const MethodResultBisection& R) {
    std::ofstream out(path);
    ensure_open(out);
    out << std::fixed << std::setprecision(12);
    out << "n,a,b,c,f(c)\n";
    for (const auto& r : R.rows) {
        out << r.n << "," << r.a << "," << r.b << "," << r.c << "," << r.fc << "\n";
    }
}

// CSV для простых итераций/Ньютона
void save_iterations_csv(const std::string& path, const MethodResultIter& R) {
    std::ofstream out(path);
    ensure_open(out);
    out << std::fixed << std::setprecision(12);
    out << "n,x,f(x),delta,residual\n";
    for (const auto& r : R.rows) {
        out << r.n << "," << r.x << "," << r.fx << "," << r.delta << "," << r.residual << "\n";
    }
}

/***********************
 * (NEW) Экспорт summary CSV
 ***********************/
void save_summary_csv(const std::string& path, const std::vector<SummaryRow>& rows) {
    std::ofstream out(path);
    ensure_open(out);
    out << std::fixed << std::setprecision(12);

    // Заголовок
    out << "method,eps,root,f(root),iterations\n";

    // Данные
    for (const auto& r : rows) {
        out << r.method << ","
            << r.eps    << ","
            << r.root   << ","
            << r.fval   << ","
            << r.iters  << "\n";
    }
}