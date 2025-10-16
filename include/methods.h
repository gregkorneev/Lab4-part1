#pragma once
#include <vector>
#include <string>
#include <functional>

struct BisectionRow {
    int    n;
    double a;
    double b;
    double c;
    double fc;
};

struct IterRow {
    int    n;
    double x;
    double fx;
    double delta;     // |x_{n+1} - x_n|
    double residual;  // |f(x_{n+1})|
};

struct MethodResultBisection {
    double root;
    int    iterations;
    double f_at_root;
    double eps;
    std::vector<BisectionRow> rows;
};

struct MethodResultIter {
    double root;
    int    iterations;
    double f_at_root;
    double eps;
    std::vector<IterRow> rows;
};

/// f(x) = x^3 - 2
double f(double x);

/// f'(x) = 3 x^2
double df(double x);

/// Метод половинного деления (бисекции)
MethodResultBisection solve_bisection(double a0, double b0, double eps, int Nmax = 200);

/// Метод простых итераций x_{n+1} = phi(x_n),
/// используем phi(x) = (x + 2/x^2)/2  (контракт в окрестности корня)
MethodResultIter solve_fixed_point(double x0, double eps, int Nmax = 200);

/// Метод Ньютона
MethodResultIter solve_newton(double x0, double eps, int Nmax = 200);

/// Экспорт таблиц в CSV
void save_bisection_csv(const std::string& path, const MethodResultBisection& R);
void save_iterations_csv(const std::string& path, const MethodResultIter& R);