#include <cmath>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>

std::deque<std::pair<double, double>>
gen_init_points(std::uint32_t points, std::function<double(double)> fn) {
    double l = -std::numbers::pi / 2.0;
    double r = std::numbers::pi / 2.0;
    // Pairs of x and f(x)
    std::deque<std::pair<double, double>> points_q;
    double step = (r - l) / (double)(points - 1);
    double cur = l;
    for (std::uint32_t i = 0; i < points; i++) {
        points_q.push_back({cur, fn(cur)});
        cur += step;
    }
    return points_q;
}
int EVALS = 0;
double f(double x) {
    double sinx = std::sin(x);
    EVALS += 1;
    return std::exp(sinx * sinx);
}

double integrate_tapezoid(std::function<double(double)> fn) {
    // Pairs of x and f(x)
    std::deque<std::pair<double, double>> iter_q = gen_init_points(2, fn);

    double I_n = -10000;
    double I_n1 = 0;
    std::cout << "=================\n";
    std::cout << "Trapezoid \n";
    std::cout << "=================\n";

    while (!(std::abs(I_n - I_n1) < 1e-6)) {
        std::deque<std::pair<double, double>> temp_q;
        I_n = I_n1;
        I_n1 = 0;
        // calculate I_n1 from simpson method
        while (iter_q.size() >= 2) {
            auto p1 = iter_q.front();
            iter_q.pop_front();
            auto p2 = iter_q.front();
            // We dont remove last node becouse next interval uses it

            I_n1 += ((p2.first - p1.first) / 2.0) * (p1.second + p2.second);

            // Add new intervals to second stack
            temp_q.push_back(p1);
            auto m1 = (p1.first + p2.first) / 2.0;
            temp_q.push_back({m1, fn(m1)});
            // We dont add r element here as it will be added with next range of
            // elements
        }
        // Add last node of range to maintain correct number of points
        temp_q.push_back(iter_q.front());
        // swap deque
        iter_q = std::move(temp_q);
    }
    std::cout << "Integrated function evals " << iter_q.size() << "\n";
    return I_n1;
}

double integrate_simpson(std::function<double(double)> fn) {
    // Pairs of x and f(x)
    std::deque<std::pair<double, double>> iter_q = gen_init_points(3, fn);

    double I_n = -10000;
    double I_n1 = 0;
    std::cout << "=================\n";
    std::cout << "Simpson \n";
    std::cout << "=================\n";

    while (!(std::abs(I_n - I_n1) < 1e-6)) {
        std::deque<std::pair<double, double>> temp_q;
        I_n = I_n1;
        I_n1 = 0;
        // calculate I_n1 from simpson method
        while (iter_q.size() >= 2) {
            auto p1 = iter_q.front();
            iter_q.pop_front();
            auto p2 = iter_q.front();
            iter_q.pop_front();
            auto p3 = iter_q.front();
            // We dont remove last node becouse next interval uses it

            I_n1 += ((p3.first - p1.first) / 6.0) *
                    (p1.second + 4 * p2.second + p3.second);

            // Add new intervals to second stack
            temp_q.push_back(p1);
            auto m1 = (p1.first + p2.first) / 2.0;
            temp_q.push_back({m1, fn(m1)});
            temp_q.push_back(p2);
            auto m2 = (p3.first + p2.first) / 2.0;
            temp_q.push_back({m2, fn(m2)});
            // We dont add r element here as it will be added with next range of
            // elements
        }
        // Add last node of range to maintain correct number of points
        temp_q.push_back(iter_q.front());
        // swap deque
        iter_q = std::move(temp_q);
    }
    std::cout << "Integrated function evals " << iter_q.size() << "\n";
    return I_n1;
}

double integrate_38(std::function<double(double)> fn) {
    // Pairs of x and f(x)
    std::deque<std::pair<double, double>> iter_q = gen_init_points(4, fn);

    double I_n = -10000;
    double I_n1 = 0;
    std::cout << "=================\n";
    std::cout << "Reguła 3/8 \n";
    std::cout << "=================\n";

    while (!(std::abs(I_n - I_n1) < 1e-6)) {
        std::deque<std::pair<double, double>> temp_q;
        I_n = I_n1;
        I_n1 = 0;
        // calculate I_n1 from simpson method
        while (iter_q.size() >= 4) {
            auto p1 = iter_q.front();
            iter_q.pop_front();
            auto p2 = iter_q.front();
            iter_q.pop_front();
            auto p3 = iter_q.front();
            iter_q.pop_front();
            auto p4 = iter_q.front();
            // We dont remove last node becouse next interval uses it

            I_n1 += ((p4.first - p1.first) / 8.0) *
                    (p1.second + p4.second + 3 * p2.second + 3 * p3.second);

            // Add new intervals to second stack
            temp_q.push_back(p1);
            auto m1 = (p1.first + p2.first) / 2.0;
            temp_q.push_back({m1, fn(m1)});
            temp_q.push_back(p2);
            auto m2 = (p2.first + p3.first) / 2.0;
            temp_q.push_back({m2, fn(m2)});
            temp_q.push_back(p3);
            auto m3 = (p3.first + p4.first) / 2.0;
            temp_q.push_back({m3, fn(m3)});
            // We dont add p4 element here as it will be added with next range
            // of elements
        }
        // Add last node of range to maintain correct number of points
        temp_q.push_back(iter_q.front());
        // swap deque
        iter_q = std::move(temp_q);
    }
    std::cout << "Integrated function evals " << iter_q.size() << "\n";
    return I_n1;
}

auto main() -> int {
    std::cout << std::setprecision(12) << integrate_tapezoid(f) << "\n";
    std::cout << "Function was calculated " << EVALS << " times \n";
    EVALS = 0;
    std::cout << std::setprecision(12) << integrate_simpson(f) << "\n";
    std::cout << "Function was calculated " << EVALS << " times \n";
    EVALS = 0;

    std::cout << std::setprecision(12) << integrate_38(f) << "\n";
    std::cout << "Function was calculated " << EVALS << " times \n";

    return 0;
}