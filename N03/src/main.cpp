#include "Eigen/Core"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

void save_to_file(std::string filename,
                  std::pair<std::string, std::string> columns,
                  std::vector<std::pair<std::uint32_t, double>> points) {
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

class ThreeDiagMat {
    /// ud - upper diag 1 - N-1 (a_N is ignored)
    Eigen::VectorXd ud;
    /// md - diag 1 - N
    Eigen::VectorXd md;
    /// ld - lower diag 2 - N (c_1 is ignored)
    Eigen::VectorXd ld;

    public:
    ThreeDiagMat(Eigen::VectorXd &&x, Eigen::VectorXd &&y, Eigen::VectorXd &&z)
        : ud(std::move(x)), md(std::move(y)), ld(std::move(z)) {}

    Eigen::VectorXd const operator*(Eigen::VectorXd const &rhs) {
        long N = md.size();
        Eigen::VectorXd sol(N);

        // First row, no lower diag element
        sol[0] = md[0] * rhs[0] + ud[0] * rhs[1];
        // 2-(N-1) rows, all three elements
        for (int i = 1; i < N - 1; ++i) {
            sol[i] = ld[i] * rhs[i - 1] + md[i] * rhs[i] + ud[i] * rhs[i + 1];
        }
        // Last row, no upper diag element
        sol[N - 1] = ld[N - 1] * rhs[N - 2] + md[N - 1] * rhs[N - 1];
        return sol;
    }

    // Adds c to all elements of diagonal. Equiv of mat = mat + c*I
    void add_diag(double c) {
        for (int i = 0; i < md.size(); i++) {
            md[i] += c;
        }
    }

    long size() { return md.size(); }
    friend class ThomasSolver;
};

class ThomasSolver {
    // Eigen::VectorXd x;
    Eigen::VectorXd y;
    Eigen::VectorXd z;
    Eigen::VectorXd c;
    /// Creates Thomas Solver for three diag matrix
    /// x upper diag 1 - N-1 (a_N is ignored)
    /// y diag 1 - N
    /// z lower diag 2 - N (c_1 is ignored)
    public:
    // Based on: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    ThomasSolver(Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd z)
        : y(y), z(z) {
        long N = this->y.size();
        this->c = Eigen::VectorXd(N);

        this->c(0) = x(0) / this->y(0);
        for (long i = 1; i < this->y.size(); i++) {
            this->c(i) = x(i) / (this->y(i) - this->z(i) * this->c(i - 1));
        }
    }

    ThomasSolver(ThreeDiagMat &mat) : ThomasSolver(mat.ud, mat.md, mat.ld) {}

    Eigen::VectorXd solve(const Eigen::VectorXd &b) const {
        auto N = this->y.size();
        assert(b.size() == N);
        Eigen::VectorXd d = Eigen::VectorXd::Zero(N);

        d(0) = b(0) / y(0);
        // Calculate coefficients in forwoard sweep
        for (long i = 1; i < this->y.size(); i++) {
            d(i) = (b(i) - this->z(i) * d(i - 1)) /
                   (this->y(i) - this->z(i) * this->c(i - 1));
        }

        auto x = Eigen::VectorXd(N);
        // Backsubstitute for x_n
        x(this->y.size() - 1) = d(N - 1);
        for (long i = N - 2; i >= 0; i--) {
            x(i) = d(i) - c(i) * x(i + 1);
        }

        return x;
    }
};

/// Generate matrix
ThreeDiagMat gen_mat(std::uint32_t N) {
    double h = 20 / ((double)N - 1.0);
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(N, N);

    Eigen::VectorXd ud(N);
    Eigen::VectorXd md(N);
    Eigen::VectorXd ld(N);

    for (long i = 0; i < N; i++) {
        double v = ((double)i * h - 10.0) * ((double)i * h - 10.0);
        md(i) = -2.0 / (h * h) + v;
        ud(i) = 1.0 / (h * h);
        ld(i) = 1.0 / (h * h);
    }
    return ThreeDiagMat(std::move(ud), std::move(md), std::move(ld));
}

Eigen::VectorXd proj(Eigen::VectorXd x, Eigen::VectorXd y) {
    Eigen::VectorXd proj_vec = x.dot(y) * x;
    return proj_vec;
}

// Calculate N bigest eigenvalues using power method
std::vector<std::pair<double, Eigen::VectorXd>>
calc_biggest_eigen_vals_power(ThreeDiagMat mat, std::uint32_t vals_count) {
    vals_count = std::min((std::uint32_t)mat.size(), vals_count);
    std::vector<std::pair<double, Eigen::VectorXd>> solutions{};

    double eps = 1e-8;
    std::uint32_t max_iter = 5000;

    // Loop for every eigenvalue
    for (std::uint32_t i = 0; i < vals_count; i++) {
        Eigen::VectorXd y_n = Eigen::VectorXd::Random(mat.size());
        Eigen::VectorXd z_n = Eigen::VectorXd(mat.size());
        y_n.normalize();

        std::uint32_t iter = 0;
        std::vector<std::pair<std::uint32_t, double>> approx{};

        do {
            z_n = mat * y_n;
            for (auto &[val, vec] : solutions) {
                z_n = z_n - proj(vec, z_n);
            }

            approx.push_back({iter, y_n.dot(mat * y_n)});
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

        // Make approximations dependant on final eigenvalue
        for (auto &[_, a] : approx) {
            a = std::abs(std::abs(a) - std::abs(val));
        }
        save_to_file("./data/biggest_iter_power" + std::to_string(i) + ".csv",
                     {"N", "Approx"}, approx);
    }

    return solutions;
}

// Calculate N smallest eigenvalues using power method
std::vector<std::pair<double, Eigen::VectorXd>>
calc_smallest_eigen_vals_power(ThreeDiagMat mat, std::uint32_t vals_count) {
    vals_count = std::min((std::uint32_t)mat.size(), vals_count);
    std::vector<std::pair<double, Eigen::VectorXd>> solutions{};

    double eps = 1e-8;

    std::uint32_t max_iter = 5000;
    ThomasSolver thomas(mat);

    // Loop for every eigenvalue
    for (std::uint32_t i = 0; i < vals_count; i++) {
        Eigen::VectorXd y_n = Eigen::VectorXd::Random(mat.size());
        Eigen::VectorXd z_n = Eigen::VectorXd(mat.size());
        y_n.normalize();

        std::uint32_t iter = 0;
        std::vector<std::pair<std::uint32_t, double>> approx{};

        do {
            z_n = thomas.solve(y_n);
            for (auto &[val, vec] : solutions) {
                z_n = z_n - proj(vec, z_n);
            }

            // approx.push_back({iter, 1.0 / (y_n.dot(thomas.solve(y_n)))});
            approx.push_back({iter, 1.0 / z_n.norm()});
            z_n.normalize();

            // It make sense
            if (std::abs(y_n.dot(z_n)) > 1.0 - eps)
                break;

            y_n = z_n;
            iter++;
        } while (iter < max_iter);
        // Calculate eigen_val and push vector
        double val = 1.0 / (y_n.dot(thomas.solve(y_n)));
        solutions.push_back(std::make_pair(val, y_n));

        // Make approximations dependant on final eigenvalue
        for (auto &[_, a] : approx) {
            a = std::abs(std::abs(a) - std::abs(val));
        }
        save_to_file("./data/samllest_iter_power" + std::to_string(i) + ".csv",
                     {"N", "Approx"}, approx);
    }

    return solutions;
}

// Calculate N smallest eigenvalues using rayleiga method
std::vector<std::pair<double, Eigen::VectorXd>>
calc_smallest_eigen_vals_ray(ThreeDiagMat mat, std::uint32_t vals_count) {
    vals_count = std::min((std::uint32_t)mat.size(), vals_count);
    std::vector<std::pair<double, Eigen::VectorXd>> solutions{};

    double eps = 1e-8;
    std::uint32_t max_iter = 1000;

    // Loop for every eigenvalue
    for (std::uint32_t i = 0; i < vals_count; i++) {
        Eigen::VectorXd y_n = Eigen::VectorXd::Random(mat.size());
        Eigen::VectorXd z_n = Eigen::VectorXd(mat.size());
        y_n.normalize();

        std::uint32_t iter = 0;
        std::vector<std::pair<std::uint32_t, double>> approx{};

        double sigma = 0;
        do {
            if (iter < 12) {
                ThomasSolver thomas(mat);
                z_n = thomas.solve(y_n);
            } else {
                sigma = y_n.dot(mat * y_n);

                ThreeDiagMat shifted = mat;
                shifted.add_diag(-sigma);
                ThomasSolver thomas(shifted);
                z_n = thomas.solve(y_n);
            }

            for (auto &[val, vec] : solutions) {
                z_n = z_n - proj(vec, z_n);
            }

            // approx.push_back({iter, sigma});
            approx.push_back({iter, 1.0 / z_n.norm()});
            z_n.normalize();

            if (std::abs(y_n.dot(z_n)) > 1.0 - eps)
                break;

            y_n = z_n;
            iter++;
        } while (iter < max_iter);

        double val = y_n.dot(mat * y_n);
        solutions.push_back(std::make_pair(val, y_n));

        // Make approximations dependant on final eigenvalue
        for (auto &[_, a] : approx) {
            a = std::abs(std::abs(a) - std::abs(val));
        }
        save_to_file("./data/samllest_iter_ray" + std::to_string(i) + ".csv",
                     {"N", "Approx"}, approx);
    }

    return solutions;
}

// Calculate N bigest eigenvalues using rayleiga method
std::vector<std::pair<double, Eigen::VectorXd>>
calc_biggest_eigen_vals_ray(ThreeDiagMat mat, std::uint32_t vals_count) {
    vals_count = std::min((std::uint32_t)mat.size(), vals_count);
    std::vector<std::pair<double, Eigen::VectorXd>> solutions{};

    double eps = 1e-8;
    std::uint32_t max_iter = 4000;

    // Loop for every eigenvalue
    for (std::uint32_t i = 0; i < vals_count; i++) {
        Eigen::VectorXd y_n = Eigen::VectorXd::Random(mat.size());
        Eigen::VectorXd z_n = Eigen::VectorXd(mat.size());
        y_n.normalize();
        std::uint32_t iter = 0;
        std::vector<std::pair<std::uint32_t, double>> approx{};

        double sigma = 0;
        do {
            if (iter < 16) {
                z_n = mat * y_n;
                sigma = y_n.dot(mat * y_n);
            } else {
                sigma = y_n.dot(mat * y_n);
                ThreeDiagMat shifted = mat;
                double perturbation = (iter % 2 == 0) ? 1e-10 : -1e-10;
                shifted.add_diag(-(sigma + perturbation) * 1.01);
                ThomasSolver thomas(shifted);
                z_n = thomas.solve(y_n);
            }

            for (auto &[val, vec] : solutions) {
                z_n -= vec * (vec.dot(z_n));
            }
            z_n.normalize();

            if (std::abs(y_n.dot(z_n)) > 1.0 - eps)
                break;

            approx.push_back({iter, sigma});
            y_n = z_n;
            iter++;
        } while (iter < max_iter);

        double val = y_n.dot(mat * y_n);
        solutions.push_back(std::make_pair(val, y_n));

        // Make approximations dependant on final eigenvalue
        for (auto &[_, a] : approx) {
            a = std::abs(std::abs(a) - std::abs(val));
        }
        save_to_file("./data/biggest_iter_ray" + std::to_string(i) + ".csv",
                     {"N", "Approx"}, approx);
    }

    return solutions;
}

auto main() -> int {
    calc_biggest_eigen_vals_power(gen_mat(500), 2);
    calc_biggest_eigen_vals_ray(gen_mat(500), 2);
    // auto ans = calc_smallest_eigen_vals_ray(gen_mat(500), 4);
    auto ans = calc_smallest_eigen_vals_power(gen_mat(500), 4);
    for (auto &[eval, evec] : ans) {
        std::cout << eval << "\n";
    }
}