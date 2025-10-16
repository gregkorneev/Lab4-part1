#pragma once

#include <vector>
#include <string>

/*
 * Заголовок объявляет:
 *  - структуры для протоколирования итераций (построчные логи);
 *  - контейнеры результатов по каждому методу;
 *  - сигнатуры функций: f(x), f'(x), три численных метода;
 *  - экспорт таблиц итераций в CSV;
 *  - (NEW) summary-структуру и экспорт общей сводной таблицы.
 *
 * Задача: решить f(x) = x^3 - 2 = 0 (вариант №12) тремя методами:
 *  1) Бисекция
 *  2) Простые итерации
 *  3) Ньютон
 * и сохранить итерации + сводку.
 */

/*** Структуры для логирования итераций ***/

// Строка протокола для метода бисекции
struct BisectionRow {
    int    n;   // номер итерации
    double a;   // левая граница интервала
    double b;   // правая граница интервала
    double c;   // середина интервала
    double fc;  // f(c)
};

// Строка протокола для методов "в одной точке" (простые итерации, Ньютон)
struct IterRow {
    int    n;         // номер итерации
    double x;         // текущее приближение (после шага)
    double fx;        // f(x) на этом шаге
    double delta;     // |x_{n+1} - x_n|
    double residual;  // |f(x_{n+1})|
};

/*** Итоговые ко``нтейнеры результатов ***/

// Результат бисекции
struct MethodResultBisection {
    double root;       // найденный корень
    int    iterations; // кол-во итераций
    double f_at_root;  // f(root)
    double eps;        // точность
    std::vector<BisectionRow> rows; // построчный лог
};

// Результат "точечных" методов (простые итерации / Ньютон)
struct MethodResultIter {
    double root;
    int    iterations;
    double f_at_root;
    double eps;
    std::vector<IterRow> rows;
};

/*** Исходная функция и производная ***/

// f(x) = x^3 - 2 (вариант №12)
double f(double x);

// f'(x) = 3x^2 (для Ньютона)
double df(double x);

/*** Численные методы ***/

// 1) Бисекция
MethodResultBisection solve_bisection(double a0, double b0, double eps, int Nmax = 200);

// 2) Простые итерации: x_{n+1} = phi(x_n), phi(x) = 0.5 * (x + 2/x^2)
MethodResultIter solve_fixed_point(double x0, double eps, int Nmax = 200);

// 3) Ньютон: x_{n+1} = x_n - f(x_n)/f'(x_n)
MethodResultIter solve_newton(double x0, double eps, int Nmax = 200);

/*** Экспорт таблиц итераций в CSV ***/

// CSV для бисекции
void save_bisection_csv(const std::string& path, const MethodResultBisection& R);

// CSV для простых итераций/Ньютона
void save_iterations_csv(const std::string& path, const MethodResultIter& R);

/*** (NEW) Summary-таблица ***/

// Одна строка сводной таблицы по методу и точности
struct SummaryRow {
    std::string method;  // "Bisection" | "FixedPoint" | "Newton"
    double      eps;     // точность
    double      root;    // найденный корень
    double      fval;    // f(root)
    int         iters;   // кол-во итераций
};

// Экспорт сводной таблицы в CSV
void save_summary_csv(const std::string& path, const std::vector<SummaryRow>& rows);