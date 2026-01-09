#include <Eigen/LU>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

Eigen::MatrixXd gen_mat(std::uint32_t N);
Eigen::VectorXd proj(Eigen::VectorXd x, Eigen::VectorXd y);

std::vector<std::pair<double, Eigen::VectorXd>>
calc_biggest_eigen_vals_power(Eigen::MatrixXd mat, std::uint32_t vals_count);
std::vector<std::pair<double, Eigen::VectorXd>>
calc_smallest_eigen_vals_power(Eigen::MatrixXd mat, std::uint32_t vals_count);

std::vector<std::pair<double, Eigen::VectorXd>>
calc_biggest_eigen_vals_ray(Eigen::MatrixXd mat, std::uint32_t vals_count);
std::vector<std::pair<double, Eigen::VectorXd>>
calc_smallest_eigen_vals_ray(Eigen::MatrixXd mat, std::uint32_t vals_count);

auto main() -> int {
    Eigen::Vector2d v(2, 2);
    std::cout << "Hello Eigen \n" << std::endl;
    std::cout << gen_mat(5) << std::endl;
    std::cout << "==============\n";
    auto ans = calc_smallest_eigen_vals_power(gen_mat(5), 5);
    std::cout << "Vec has " << ans.size() << " eigen pairs" << std::endl;
    for (auto &[val, vec] : ans) {
        std::cout << "Eigen vector for val " << val << ":\n";
        // std::cout << vec << std::endl;
    }
}

// TODO check with sparse matrix, could be faster
Eigen::MatrixXd gen_mat(std::uint32_t N) {
    double h = 20 / ((double)N - 1.0);
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(N, N);

    for (long i = 0; i < N; i++) {
        double v = ((double)i * h - 10.0) * ((double)i * h - 10.0);
        mat(i, i) = -2.0 / (h * h) + v;
        // Above diagonal element
        if (i > 0) {
            mat(i, i - 1) = 1.0 / (h * h);
        }
        if (i < N - 1) {
            mat(i, i + 1) = 1.0 / (h * h);
        }
    }
    return mat;
}

Eigen::VectorXd proj(Eigen::VectorXd x, Eigen::VectorXd y) {
    Eigen::VectorXd proj_vec = (x.dot(y) / x.dot(x)) * x;
    return proj_vec;
}

std::vector<std::pair<double, Eigen::VectorXd>>
calc_biggest_eigen_vals_power(Eigen::MatrixXd mat, std::uint32_t vals_count) {
    vals_count = std::min((std::uint32_t)mat.rows(), vals_count);
    std::vector<std::pair<double, Eigen::VectorXd>> solutions{};
    std::cout << "Values to calculate " << vals_count << std::endl;

    double eps = 1e-8;
    std::uint32_t max_iter = 1000;

    // Loop for every eigenvalue
    for (std::uint32_t i = 0; i < vals_count; i++) {
        Eigen::VectorXd y_n = Eigen::VectorXd::Random(mat.rows());
        Eigen::VectorXd z_n = Eigen::VectorXd(mat.rows());
        y_n.normalize();
        std::uint32_t iter = 0;

        do {
            z_n = mat * y_n;
            for (auto &[val, vec] : solutions) {
                z_n = z_n - proj(vec, z_n);
            }
            z_n.normalize();

            // It make sense
            if (std::abs(y_n.dot(z_n)) > 1.0 - eps)
                break;

            y_n = z_n;
            iter++;
        } while (iter < max_iter);
        // Calculate eigen_val and push vector
        double val = y_n.dot(mat * y_n);
        solutions.push_back(std::make_pair(val, y_n));
    }

    return solutions;
}

std::vector<std::pair<double, Eigen::VectorXd>>
calc_smallest_eigen_vals_power(Eigen::MatrixXd mat, std::uint32_t vals_count) {
    vals_count = std::min((std::uint32_t)mat.rows(), vals_count);
    std::vector<std::pair<double, Eigen::VectorXd>> solutions{};
    std::cout << "Values to calculate " << vals_count << std::endl;

    double eps = 1e-8;
    std::uint32_t max_iter = 5000;
    Eigen::PartialPivLU<Eigen::MatrixXd> lu(mat);

    // Loop for every eigenvalue
    for (std::uint32_t i = 0; i < vals_count; i++) {
        Eigen::VectorXd y_n = Eigen::VectorXd::Random(mat.rows());
        Eigen::VectorXd z_n = Eigen::VectorXd(mat.rows());
        y_n.normalize();
        std::uint32_t iter = 0;

        do {
            z_n = lu.solve(y_n);
            for (auto &[val, vec] : solutions) {
                z_n = z_n - proj(vec, z_n);
            }
            z_n.normalize();

            // It make sense
            if (std::abs(y_n.dot(z_n)) > 1.0 - eps)
                break;

            y_n = z_n;
            iter++;
        } while (iter < max_iter);
        // Calculate eigen_val and push vector
        double val = 1.0 / (y_n.dot(lu.solve(y_n)));
        solutions.push_back(std::make_pair(val, y_n));
    }

    return solutions;
}
