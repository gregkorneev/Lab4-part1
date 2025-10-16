#include "methods.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

double f(double x) {
    return x*x*x - 2.0;
}

double df(double x) {
    return 3.0 * x * x;
}

MethodResultBisection solve_bisection(double a, double b, double eps, int Nmax) {
    if (a >= b) throw std::invalid_argument("bisection: a >= b");
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0.0) throw std::invalid_argument("bisection: f(a) and f(b) must have opposite signs");

    MethodResultBisection res{};
    res.eps = eps;

    for (int n = 0; n < Nmax; ++n) {
        double c = 0.5 * (a + b);
        double fc = f(c);

        res.rows.push_back({n, a, b, c, fc});

        if (std::fabs(fc) < eps || 0.5 * (b - a) < eps) {
            res.root = c;
            res.f_at_root = fc;
            res.iterations = n + 1;
            return res;
        }

        // Выбор подынтервала
        if (fa * fc < 0.0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }

    // Если достигли Nmax без критерия
    double c = 0.5 * (a + b);
    res.root = c;
    res.f_at_root = f(c);
    res.iterations = Nmax;
    return res;
}

MethodResultIter solve_fixed_point(double x0, double eps, int Nmax) {
    auto phi = [](double x) {
        // phi(x) = (x + 2/x^2)/2  — сходимость: phi'(root) = -0.5
        return 0.5 * (x + 2.0 / (x * x));
    };

    MethodResultIter res{};
    res.eps = eps;

    double x = x0;
    for (int n = 0; n < Nmax; ++n) {
        // защита от деления на 0
        if (std::fabs(x) < 1e-14) x = 1e-6;

        double xn1 = phi(x);
        double delta = std::fabs(xn1 - x);
        double r = std::fabs(f(xn1));

        res.rows.push_back({n, xn1, f(xn1), delta, r});

        if (delta < eps || r < eps) {
            res.root = xn1;
            res.f_at_root = f(xn1);
            res.iterations = n + 1;
            return res;
        }
        x = xn1;
    }

    res.root = x;
    res.f_at_root = f(x);
    res.iterations = Nmax;
    return res;
}

MethodResultIter solve_newton(double x0, double eps, int Nmax) {
    MethodResultIter res{};
    res.eps = eps;

    double x = x0;
    for (int n = 0; n < Nmax; ++n) {
        double y = f(x);
        double dy = df(x);

        if (std::fabs(dy) < 1e-14) {
            // производная близка к нулю — остановим, чтобы не взорваться
            break;
        }

        double xn1 = x - y / dy;
        double delta = std::fabs(xn1 - x);
        double r = std::fabs(f(xn1));

        res.rows.push_back({n, xn1, f(xn1), delta, r});

        if (delta < eps || r < eps) {
            res.root = xn1;
            res.f_at_root = f(xn1);
            res.iterations = n + 1;
            return res;
        }
        x = xn1;
    }

    res.root = x;
    res.f_at_root = f(x);
    res.iterations = (int)res.rows.size();
    return res;
}

static void ensure_csv_dir(std::ofstream& out) {
    if (!out.good())
        throw std::runtime_error("Cannot open CSV file for writing");
}

void save_bisection_csv(const std::string& path, const MethodResultBisection& R) {
    std::ofstream out(path);
    ensure_csv_dir(out);
    out << std::fixed << std::setprecision(12);
    out << "n,a,b,c,f(c)\n";
    for (const auto& r : R.rows) {
        out << r.n << "," << r.a << "," << r.b << "," << r.c << "," << r.fc << "\n";
    }
}

void save_iterations_csv(const std::string& path, const MethodResultIter& R) {
    std::ofstream out(path);
    ensure_csv_dir(out);
    out << std::fixed << std::setprecision(12);
    out << "n,x,f(x),delta,residual\n";
    for (const auto& r : R.rows) {
        out << r.n << "," << r.x << "," << r.fx << "," << r.delta << "," << r.residual << "\n";
    }
}