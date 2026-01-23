#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

std::vector<std::pair<double, double>>
gen_init_points(std::uint32_t points, std::function<double(double)> fn) {
    double l = -std::numbers::pi / 2.0;
    double r = std::numbers::pi / 2.0;
    // Pairs of x and f(x)
    std::vector<std::pair<double, double>> points_q;
    double step = (r - l) / (double)(points - 1);
    double cur = l;
    for (std::uint32_t i = 0; i < points; i++) {
        points_q.push_back({cur, fn(cur)});
        cur += step;
    }
    return points_q;
}
// CHECK HOW MANY TIMES FUNCTION WAS CALLED
int EVALS = 0;
double f(double x) {
    double sinx = std::sin(x);
    EVALS += 1;
    return std::exp(sinx * sinx);
}

std::vector<std::pair<double, double>>
densify_points(std::vector<std::pair<double, double>>& iter_q,
               std::function<double(double)> fn) {
    std::vector<std::pair<double, double>> temp_q;
    for (unsigned int i = 0; i < iter_q.size() - 1; i++) {
        auto p1 = iter_q[i];
        auto p2 = iter_q[i + 1];
        temp_q.push_back(p1);
        auto m1 = (p1.first + p2.first) / 2.0;
        temp_q.push_back({m1, fn(m1)});
    }
    temp_q.push_back(iter_q.back()); // add last element

    return temp_q;
}

double integrate_tapezoid(std::function<double(double)> fn) {
    // Pairs of x and f(x)
    std::vector<std::pair<double, double>> iter_q = gen_init_points(2, fn);

    double I_n = -10000;
    double I_n1 = 0;
    std::cout << "=================\n";
    std::cout << "Trapezoid \n";
    std::cout << "=================\n";

    while (!(std::abs(I_n - I_n1) < 1e-6)) {
        // Make set denser
        // We calculate new points after check if our prev approximation was
        // sufficient to remove redundant calculations of points if the next
        // iteration wouldn't happen. It technicly starting point count could be
        // sufficient but it isnt likely
        iter_q = densify_points(iter_q, fn);

        I_n = I_n1;
        I_n1 = 0;
        // calculate I_n1 from simpson method
        for (unsigned int i = 0; i < iter_q.size() - 1; i += 1) {
            auto p1 = iter_q[i];
            auto p2 = iter_q[i + 1];

            I_n1 += ((p2.first - p1.first) / 2.0) * (p1.second + p2.second);
        }
    }
    std::cout << "Integrated function evals " << iter_q.size() << "\n";

    // Open file for data write
    // Save iter_q to file to make plot
    std::fstream file;
    file.open("./data/trapezoid_points.csv", std::ios::out);
    if (!file) {
        std::cerr << "Cant open file to save data: "
                  << "./data/trapezoid_points.csv" << std::endl;
        exit(1);
    }
    file << "x" << ";" << "y";
    file << std::endl;

    for (auto& [x, y] : iter_q) {
        file << x << ";" << y << "\n";
    }
    file.close();

    return I_n1;
}

double integrate_simpson(std::function<double(double)> fn) {
    // Pairs of x and f(x)
    std::vector<std::pair<double, double>> iter_q = gen_init_points(3, fn);
    double I_n = -10000;
    double I_n1 = 0;
    std::cout << "=================\n";
    std::cout << "Simpson \n";
    std::cout << "=================\n";

    while (!(std::abs(I_n - I_n1) < 1e-6)) {
        // Make set denser
        // We calculate new points after check if our prev approximation was
        // sufficient to remove redundant calculations of points if the next
        // iteration wouldn't happen. It technicly starting point count could be
        // sufficient but it isnt likely
        iter_q = densify_points(iter_q, fn);

        I_n = I_n1;
        I_n1 = 0;
        // calculate I_n1 from simpson method
        for (unsigned int i = 0; i < iter_q.size() - 2; i += 2) {
            auto p1 = iter_q[i];
            auto p2 = iter_q[i + 1];
            auto p3 = iter_q[i + 2];

            I_n1 += ((p3.first - p1.first) / 6.0) *
                    (p1.second + 4 * p2.second + p3.second);
        }
    }
    std::cout << "Integrated function evals " << iter_q.size() << "\n";
    return I_n1;
}

double integrate_38(std::function<double(double)> fn) {
    // Pairs of x and f(x)
    std::vector<std::pair<double, double>> iter_q = gen_init_points(4, fn);

    double I_n = -10000;
    double I_n1 = 0;
    std::cout << "=================\n";
    std::cout << "Reguła 3/8 \n";
    std::cout << "=================\n";

    while (!(std::abs(I_n - I_n1) < 1e-6)) {
        // Make set denser
        // We calculate new points after check if our prev approximation was
        // sufficient to remove redundant calculations of points if the next
        // iteration wouldn't happen. It technicly starting point count could be
        // sufficient but it isnt likely
        iter_q = densify_points(iter_q, fn);

        I_n = I_n1;
        I_n1 = 0;
        // calculate I_n1 from simpson method
        for (unsigned int i = 0; i < iter_q.size() - 3; i += 3) {
            auto p1 = iter_q[i];
            auto p2 = iter_q[i + 1];
            auto p3 = iter_q[i + 2];
            auto p4 = iter_q[i + 3];

            I_n1 += ((p4.first - p1.first) / 8.0) *
                    (p1.second + p4.second + 3 * p2.second + 3 * p3.second);
        }
    }
    std::cout << "Integrated function evals " << iter_q.size() << "\n";
    return I_n1;
}

auto main() -> int {
    std::cout << std::setprecision(16) << integrate_tapezoid(f) << "\n";
    std::cout << "Function was calculated " << EVALS << " times \n";
    EVALS = 0;
    std::cout << std::setprecision(16) << integrate_simpson(f) << "\n";
    std::cout << "Function was calculated " << EVALS << " times \n";
    EVALS = 0;

    std::cout << std::setprecision(16) << integrate_38(f) << "\n";
    std::cout << "Function was calculated " << EVALS << " times \n";

    return 0;
}