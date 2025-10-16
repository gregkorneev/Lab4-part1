#include "methods.h"
#include <iostream>
#include <iomanip>

int main() {
    // Точности по заданию
    double eps1 = 1e-6;
    double eps2 = 1e-11;

    // Начальные данные для нашего уравнения f(x) = x^3 - 2
    // Бисекция требует интервал, где знаки разные: f(1) < 0, f(2) > 0
    double a0 = 1.0;
    double b0 = 2.0;

    // Стартовые приближения для итераций и Ньютона (просто разумные числа около корня ~1.26)
    double x0_fixed = 1.3;
    double x0_newton = 1.5;

    // Красивый вывод чисел
    std::cout << std::fixed << std::setprecision(12);

    // ====== eps = 1e-6 ======
    MethodResultBisection b1 = solve_bisection(a0, b0, eps1);
    MethodResultIter      i1 = solve_fixed_point(x0_fixed, eps1);
    MethodResultIter      n1 = solve_newton(x0_newton, eps1);

    // Сохраняем таблицы в CSV (они будут созданы рядом с программой в папке results, если она есть)
    // Если папки нет — просто закоммить папку results заранее, как мы делали при структуре проекта.
    save_bisection_csv("results/bisection_eps1e-6.csv", b1);
    save_iterations_csv("results/iter_eps1e-6.csv", i1);
    save_iterations_csv("results/newton_eps1e-6.csv", n1);

    // Печатаем краткие итоги для eps = 1e-6
    std::cout << "=== Results for eps = 1e-6 ===\n";
    std::cout << "Bisection: root = " << b1.root
              << ", f(root) = " << b1.f_at_root
              << ", iterations = " << b1.iterations << "\n";

    std::cout << "Fixed pt : root = " << i1.root
              << ", f(root) = " << i1.f_at_root
              << ", iterations = " << i1.iterations << "\n";

    std::cout << "Newton   : root = " << n1.root
              << ", f(root) = " << n1.f_at_root
              << ", iterations = " << n1.iterations << "\n\n";

    // ====== eps = 1e-11 ======
    MethodResultBisection b2 = solve_bisection(a0, b0, eps2);
    MethodResultIter      i2 = solve_fixed_point(x0_fixed, eps2);
    MethodResultIter      n2 = solve_newton(x0_newton, eps2);

    // Сохраняем таблицы в CSV для второй точности
    save_bisection_csv("results/bisection_eps1e-11.csv", b2);
    save_iterations_csv("results/iter_eps1e-11.csv", i2);
    save_iterations_csv("results/newton_eps1e-11.csv", n2);

    // Печатаем краткие итоги для eps = 1e-11
    std::cout << "=== Results for eps = 1e-11 ===\n";
    std::cout << "Bisection: root = " << b2.root
              << ", f(root) = " << b2.f_at_root
              << ", iterations = " << b2.iterations << "\n";

    std::cout << "Fixed pt : root = " << i2.root
              << ", f(root) = " << i2.f_at_root
              << ", iterations = " << i2.iterations << "\n";

    std::cout << "Newton   : root = " << n2.root
              << ", f(root) = " << n2.f_at_root
              << ", iterations = " << n2.iterations << "\n";

    // Подсказка, где лежат CSV
    std::cout << "\nCSV сохранены в папке results/\n";

    return 0;
}