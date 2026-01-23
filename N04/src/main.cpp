#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

/// Generates points as (x, f(x)) piars with x that  are equally spaced
std::vector<std::pair<long double, long double>>
gen_equal_points(std::uint32_t N, std::function<long double(long double)> f) {
    std::vector<std::pair<long double, long double>> points{};
    for (long double i = 0; i < (long double)N; i++) {
        long double x = -5 + 10 * i / ((long double)N - 1.0);
        points.push_back(std::make_pair(x, f(x)));
    }
    return points;
}

/// Generates points as (x, f(x)) piars with x that  are more dense at the ends
/// of the range
std::vector<std::pair<long double, long double>>
gen_czebyszew_points(std::uint32_t N,
                     std::function<long double(long double)> f) {
    std::vector<std::pair<long double, long double>> points{};
    for (long double i = 0; i < (long double)N; i++) {
        long double x = -5 * std::cos(i / (long double)(N - 1) * M_PI);
        points.push_back(std::make_pair(x, f(x)));
    }
    return points;
}

/// Interpolated function
long double fn(long double x) { return 1.0 / (1 + x * x); }

void save_to_file(std::string filename,
                  std::pair<std::string, std::string> columns,
                  std::vector<std::pair<long double, long double>> points) {
    // Open file for data write
    std::fstream file;
    file.open(filename, std::ios::out);
    if (!file) {
        std::cerr << "Cant open file to save data: " << filename << std::endl;
        exit(1);
    }
    file << columns.first << ";" << columns.second;
    file << std::endl;

    for (auto &[x, y] : points) {
        file << x << ";" << y << "\n";
    }
    file.close();
}

/// Calculates largrange coeficiens
std::vector<long double> calculate_l_coeficients(
    const std::vector<std::pair<long double, long double>> &interpolation_nodes,
    long double x) {
    std::uint32_t n = interpolation_nodes.size();
    std::vector<long double> coef{};

    for (std::uint32_t i = 0; i < n; i++) {
        long double li = 1.0;
        auto xi = interpolation_nodes[i].first;
        for (std::uint32_t j = 0; j < n; j++) {
            if (i != j) {
                long double xj = interpolation_nodes[j].first;
                li *= (x - xj) / (xi - xj);
            }
        }
        coef.push_back(li);
    }
    return coef;
}

/// calculates coeficients of interpolated polynomial
std::vector<long double> calculate_polynomial_coeficients(
    std::vector<std::pair<long double, long double>> &interpolation_nodes) {
    // If we have 0 and f(0) in our interpolation nodes
    for (std::uint32_t i = 0; i < interpolation_nodes.size(); i++) {
        if (interpolation_nodes[i].first == 0.0) {
            // we have idx of 0 node
            // we need to calculate lagrange interpolation in point != 0
            // replace 0 f(0) with new x f(x)
            long double new_x = interpolation_nodes[i + 1].first / 2;
            std::vector<long double> l =
                calculate_l_coeficients(interpolation_nodes, new_x);
            long double f_x = 0.0;
            for (std::uint32_t j = 0; j < l.size(); j++) {
                f_x += l[j] * interpolation_nodes[j].second;
            }
            interpolation_nodes[i] = std::make_pair(new_x, f_x);
            break;
        }
    }

    std::vector<long double> l =
        calculate_l_coeficients(interpolation_nodes, 0);

    std::vector<long double> fk{};
    // Populate fk for first iteration
    for (auto &[_, fx] : interpolation_nodes) {
        fk.push_back(fx);
    }

    std::vector<long double> a{};
    for (std::uint32_t i = 0; i < interpolation_nodes.size(); i++) {
        long double ak = 0;
        for (std::uint32_t j = 0; j < interpolation_nodes.size(); j++) {
            ak += l[j] * fk[j];
        }
        a.push_back(ak);
        // recalculate fk
        for (std::uint32_t j = 0; j < interpolation_nodes.size(); j++) {
            long double fj = (fk[j] - ak) / interpolation_nodes[j].first;
            fk[j] = fj;
        }
    }
    return a;
}

/// Evaluate polynomial using horner algorithm
long double horner_evaluate(long double x,
                            const std::vector<long double> &coef) {
    long double P = coef.back();
    for (std::uint32_t i = 1; i < coef.size(); i++) {
        P = (P * x) + coef[coef.size() - 1 - i];
    }
    return P;
}

std::pair<long, long double> find_min_error(
    std::function<std::vector<std::pair<long double, long double>>(
        std::uint32_t, std::function<long double(long double)>)>
        node_gen_fn,

    std::function<std::vector<std::pair<long double, long double>>(
        std::vector<std::pair<long double, long double>> &, std::uint32_t)>
        point_gen_fn,

    std::string error_file = "") {
    // gen dense points form normal function
    auto points = gen_equal_points(1024, fn);
    std::vector<double> errors{};
    for (std::uint32_t n = 2; n < 256; n += 1) {
        auto inter_nodes = node_gen_fn(n, fn);

        auto inter_points = point_gen_fn(inter_nodes, 1024);

        long double max_abs_diff = -1.0;
        for (std::uint32_t i = 0; i < inter_points.size(); i++) {
            max_abs_diff =
                std::max(max_abs_diff,
                         std::abs(points[i].second - inter_points[i].second));
        }
        std::cout << "N=" << n << " err=" << max_abs_diff << "\n";
        errors.push_back(max_abs_diff);
    }
    // Save errors to file
    if (error_file != "") {
        // TODO save to file
        std::fstream file;
        file.open(error_file, std::ios::out);
        if (!file) {
            std::cerr << "Cant open file to save data: " << error_file
                      << std::endl;
            exit(1);
        }
        file << "N" << ";" << "Error";
        file << std::endl;

        std::uint32_t N = 2;
        for (auto &err : errors) {
            file << N << ";" << err << "\n";
            N++;
        }
        file.close();
    }

    // find minimal error value and return it with its index
    auto min_it = std::min_element(errors.begin(), errors.end());
    return std::make_pair(std::distance(errors.begin(), min_it) + 2, *min_it);
}

void gen_optimal_polynomial_interpolation(
    std::function<std::vector<std::pair<long double, long double>>(
        std::uint32_t, std::function<long double(long double)>)>
        gen_point_fn,
    std::string filename) {
    auto dense_polynomial_point_fn =
        [](std::vector<std::pair<long double, long double>> &inter_nodes,
           std::uint32_t point_count) {
            auto coefs = calculate_polynomial_coeficients(inter_nodes);
            auto fun = [&coefs](long double x) {
                return horner_evaluate(x, coefs);
            };
            auto inter_points = gen_equal_points(point_count, fun);
            return inter_points;
        };

    auto min_err = find_min_error(gen_point_fn, dense_polynomial_point_fn,
                                  "./data/errors_" + filename);

    auto inter_points = gen_point_fn(min_err.first, fn);
    auto coefs = calculate_polynomial_coeficients(inter_points);
    auto fun = [&coefs](long double x) { return horner_evaluate(x, coefs); };

    save_to_file("./data/nodes_" + filename, {"x", "y"}, inter_points);
    save_to_file("./data/" + filename, {"x", "y"}, gen_equal_points(256, fun));
}

std::vector<std::pair<long double, long double>>
direct_lagrange(std::vector<std::pair<long double, long double>> inter_nodes,
                std::uint32_t dense_point_count) {
    // auto inter_nodes = gen_point_fn(inter_nodes_count, fn);
    std::uint32_t inter_nodes_count = inter_nodes.size();

    std::vector<long double> denominators{};
    for (std::uint32_t i = 0; i < inter_nodes_count; i++) {
        long double denominator = 1.0;
        for (std::uint32_t j = 0; j < inter_nodes_count; j++) {
            if (i != j) {
                denominator *= inter_nodes[i].first - inter_nodes[j].first;
            }
        }
        denominators.push_back(denominator);
    }

    auto f_approx = [&denominators, &inter_nodes](long double x) {
        long double omega = 1.0;
        for (auto &[x_i, y_i] : inter_nodes) {
            omega *= x - x_i;
            if (std::abs(x - x_i) < 1e-18L) { // epsilon dla long double
                return y_i;
            }
        }

        long double y = 0.0;
        for (std::uint32_t i = 0; i < inter_nodes.size(); i++) {
            y += (omega / (x - inter_nodes[i].first)) / denominators[i] *
                 inter_nodes[i].second;
        }

        return y;
    };

    auto dense_points = gen_equal_points(dense_point_count, f_approx);

    return dense_points;
}

void gen_optimal_direct_interpolation(
    std::function<std::vector<std::pair<long double, long double>>(
        std::uint32_t, std::function<long double(long double)>)>
        gen_point_fn,
    std::string filename) {
    auto dense_direct_point_fn =
        [](std::vector<std::pair<long double, long double>> &inter_nodes,
           std::uint32_t point_count) {
            return direct_lagrange(inter_nodes, point_count);
        };

    auto min_err = find_min_error(gen_point_fn, dense_direct_point_fn,
                                  "./data/errors_" + filename);

    auto inter_points = gen_point_fn(min_err.first, fn);

    save_to_file("./data/nodes_" + filename, {"x", "y"}, inter_points);
    save_to_file("./data/" + filename, {"x", "y"},
                 direct_lagrange(inter_points, 512));
}

auto main() -> int {
    gen_optimal_polynomial_interpolation(gen_equal_points,
                                         "equal_polynomial.csv");
    gen_optimal_polynomial_interpolation(gen_czebyszew_points,
                                         "czybyszew_polynomial.csv");

    gen_optimal_direct_interpolation(gen_equal_points, "equal_direct.csv");
    gen_optimal_direct_interpolation(gen_czebyszew_points,
                                     "czebyszew_direct.csv");

    //  Chart for runge fenomenon
    auto inter_points = gen_equal_points(12, fn);
    auto coefs = calculate_polynomial_coeficients(inter_points);
    auto fun = [&coefs](long double x) { return horner_evaluate(x, coefs); };

    save_to_file("./data/nodes_runge_polynomial.csv", {"x", "y"}, inter_points);
    save_to_file("./data/runge_polynomial.csv", {"x", "y"},
                 gen_equal_points(256, fun));

    return 0;
}

std::pair<long, long double> find_min_errorr(
    std::function<std::vector<std::pair<long double, long double>>(
        std::uint32_t, std::function<long double(long double)>)>
        node_gen_fn,

    std::function<std::vector<std::pair<long double, long double>>(
        std::vector<std::pair<long double, long double>> &, std::uint32_t)>
        point_gen_fn,

    std::string error_file = "");