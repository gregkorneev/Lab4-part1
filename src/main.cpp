#include "methods.h"

#include <iostream>     // cout, cerr
#include <iomanip>      // setprecision, fixed
#include <filesystem>   // create_directories
#include <vector>       // std::vector
#include <string>       // std::string

/*
 * Главный сценарий лабораторной:
 *  - создаём папку results;
 *  - задаём две точности eps (1e-6 и 1e-11);
 *  - считаем тремя методами (бисекция, простые итерации, Ньютон);
 *  - сохраняем детальные таблицы итераций CSV;
 *  - (NEW) печатаем сводную ASCII-таблицу в консоль;
 *  - (NEW) сохраняем общую summary-таблицу в results/summary.csv.
 *
 * Уравнение: f(x) = x^3 - 2 = 0 (вариант №12)
 */

static void print_summary_ascii(const std::vector<SummaryRow>& S) {
    // Красивый фиксированный формат для консоли
    std::cout << "\n=== SUMMARY (by method and eps) ===\n";
    std::cout << std::left
              << std::setw(14) << "Method"
              << std::setw(12) << "eps"
              << std::setw(18) << "root"
              << std::setw(18) << "f(root)"
              << std::setw(12) << "iterations"
              << "\n";

    std::cout << std::string(14+12+18+18+12, '-') << "\n";

    for (const auto& r : S) {
        std::cout << std::left
                  << std::setw(14) << r.method
                  << std::setw(12) << std::setprecision(12) << std::fixed << r.eps
                  << std::setw(18) << std::setprecision(12) << std::fixed << r.root
                  << std::setw(18) << std::setprecision(12) << std::fixed << r.fval
                  << std::setw(12) << r.iters
                  << "\n";
    }
    std::cout << "(also saved to results/summary.csv)\n";
}

int main() {
    try {
        // 1) Гарантируем наличие папки результатов (куда положим CSV)
        std::filesystem::create_directories("results");

        // 2) Две точности по заданию
        const double eps1 = 1e-6;
        const double eps2 = 1e-11;

        // 3) Исходные параметры
        double a0 = 1.0, b0 = 2.0; // для бисекции (меняющиеся знаки: f(1)<0, f(2)>0)
        double x0_iter = 1.3;      // старт для простых итераций
        double x0_newt = 1.5;      // старт для Ньютона

        // Формат чисел в консоли
        std::cout << std::fixed << std::setprecision(12);

        // ===== Серия 1: eps = 1e-6 =====
        MethodResultBisection Rb1 = solve_bisection(a0, b0, eps1);
        MethodResultIter      Ri1 = solve_fixed_point(x0_iter, eps1);
        MethodResultIter      Rn1 = solve_newton(x0_newt, eps1);

        save_bisection_csv("results/bisection_eps1e-6.csv", Rb1);
        save_iterations_csv("results/iter_eps1e-6.csv",      Ri1);
        save_iterations_csv("results/newton_eps1e-6.csv",    Rn1);

        std::cout << "=== Results for eps = 1e-6 ===\n";
        std::cout << "Bisection: root = " << Rb1.root
                  << ", f(root) = " << Rb1.f_at_root
                  << ", iterations = " << Rb1.iterations << "\n";
        std::cout << "Fixed pt : root = " << Ri1.root
                  << ", f(root) = " << Ri1.f_at_root
                  << ", iterations = " << Ri1.iterations << "\n";
        std::cout << "Newton   : root = " << Rn1.root
                  << ", f(root) = " << Rn1.f_at_root
                  << ", iterations = " << Rn1.iterations << "\n\n";

        // ===== Серия 2: eps = 1e-11 =====
        MethodResultBisection Rb2 = solve_bisection(a0, b0, eps2);
        MethodResultIter      Ri2 = solve_fixed_point(x0_iter, eps2);
        MethodResultIter      Rn2 = solve_newton(x0_newt, eps2);

        save_bisection_csv("results/bisection_eps1e-11.csv", Rb2);
        save_iterations_csv("results/iter_eps1e-11.csv",      Ri2);
        save_iterations_csv("results/newton_eps1e-11.csv",    Rn2);

        std::cout << "=== Results for eps = 1e-11 ===\n";
        std::cout << "Bisection: root = " << Rb2.root
                  << ", f(root) = " << Rb2.f_at_root
                  << ", iterations = " << Rb2.iterations << "\n";
        std::cout << "Fixed pt : root = " << Ri2.root
                  << ", f(root) = " << Ri2.f_at_root
                  << ", iterations = " << Ri2.iterations << "\n";
        std::cout << "Newton   : root = " << Rn2.root
                  << ", f(root) = " << Rn2.f_at_root
                  << ", iterations = " << Rn2.iterations << "\n";

        // 4) (NEW) Сформируем общую summary-таблицу по всем 3 методам и 2 eps
        std::vector<SummaryRow> summary;
        summary.push_back({"Bisection",  eps1, Rb1.root, Rb1.f_at_root, Rb1.iterations});
        summary.push_back({"FixedPoint", eps1, Ri1.root, Ri1.f_at_root, Ri1.iterations});
        summary.push_back({"Newton",     eps1, Rn1.root, Rn1.f_at_root, Rn1.iterations});
        summary.push_back({"Bisection",  eps2, Rb2.root, Rb2.f_at_root, Rb2.iterations});
        summary.push_back({"FixedPoint", eps2, Ri2.root, Ri2.f_at_root, Ri2.iterations});
        summary.push_back({"Newton",     eps2, Rn2.root, Rn2.f_at_root, Rn2.iterations});

        // 5) Печать ASCII-таблицы и сохранение в CSV
        print_summary_ascii(summary);
        save_summary_csv("results/summary.csv", summary);

        std::cout << "\nПодробные таблицы итераций — в ./results/ (относительно каталога запуска)\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}