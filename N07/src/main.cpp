#include "Eigen/Core"
#include <fstream>
#include <iostream>

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
    ThomasSolver(ThreeDiagMat &mat) : ThomasSolver(mat.ud, mat.md, mat.ld) {}
    // Based on: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    ThomasSolver(Eigen::VectorXd x, Eigen::VectorXd y, Eigen::VectorXd z)
        : y(y), z(z) {
        long N = this->y.size();
        this->c = Eigen::VectorXd::Zero(N);

        this->c(0) = x(0) / this->y(0);
        for (long i = 1; i < N - 1; i++) {
            this->c(i) = x(i) / (this->y(i) - this->z(i) * this->c(i - 1));
        }
    }

    Eigen::VectorXd solve(const Eigen::VectorXd &b) const {
        auto N = this->y.size();
        Eigen::VectorXd d = Eigen::VectorXd::Zero(N);

        d(0) = b(0) / y(0);
        for (long i = 1; i < N; i++) {
            d(i) = (b(i) - this->z(i) * d(i - 1)) /
                   (this->y(i) - this->z(i) * this->c(i - 1));
        }

        auto sol = Eigen::VectorXd(N);
        sol(N - 1) = d(N - 1);
        for (long i = N - 2; i >= 0; i--) {
            sol(i) = d(i) - c(i) * sol(i + 1);
        }
        return sol;
    }
};

/// Generate matrix
ThreeDiagMat gen_jacobian(Eigen::VectorXd u) {
    const std::uint32_t N = u.size();
    // In task we have N+1 functions so our N is N+1 from task.
    double h = 20.0 / ((double)N - 2.0);
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(N, N);

    Eigen::VectorXd ud(N); // uper diag
    Eigen::VectorXd md(N); // diag
    Eigen::VectorXd ld(N); // lower diag

    double h2 = h * h;
    for (long i = 0; i < N; i++) {
        md(i) = 2 / (h2) + 6.0 * u[i] * u[i] - 2.0;
        ud(i) = -1.0 / (h2);
        ld(i) = -1.0 / (h2);
    }
    // Fix first and last row
    ud[0] = 0;
    md[0] = 1;

    ld[N - 1] = 0;
    md[N - 1] = 1;

    return ThreeDiagMat(std::move(ud), std::move(md), std::move(ld));
}

Eigen::VectorXd calculate_g(const Eigen::VectorXd &u) {
    const std::uint32_t N = u.size();
    double h = 20.0 / ((double)N - 2.0);

    Eigen::VectorXd g(N);
    for (std::uint32_t i = 1; i < N - 1; i++) {
        g[i] = -(u[i - 1] - 2 * u[i] + u[i + 1]) / (h * h) +
               2 * u[i] * (u[i] * u[i] - 1.0);
    }
    g[0] = u[0];
    g[N - 1] = u[N - 1] - 1.0;
    return g;
}

Eigen::VectorXd damped_newton_solve() {
    std::uint32_t N = 999;
    Eigen::VectorXd u = Eigen::VectorXd::Zero(N);

    for (std::uint32_t i = 0; i < N; ++i) {
        u[i] = 1;
    }
    u[0] = 0;

    const double epsilon = 1e-6;
    const std::uint32_t iter_max = 500;
    std::uint32_t iter = 0;
    Eigen::VectorXd g;

    while (iter < iter_max) {
        g = calculate_g(u);
        // Check if our residual is smaller then tolerance
        double current_norm = g.norm();
        if (current_norm < epsilon)
            break;

        ThreeDiagMat J = gen_jacobian(u);
        ThomasSolver thomas(J);
        auto delta_u = thomas.solve(g);

        // u -= delta_u;

        Eigen::VectorXd u_test;

        bool accepted = false;
        double w = 1.0;
        while (w > 1e-4) {
            u_test = u - w * delta_u;

            Eigen::VectorXd g_test = calculate_g(u_test);
            double test_norm = g_test.norm();

            if (test_norm < current_norm) {
                u = u_test;
                accepted = true;
                break;
            } else {
                w /= 2.0;
            }
        }

        if (!accepted) {
            u = u_test;
        }

        iter++;
    }
    // std::cout << u;
    return u;
}

void save_to_file(std::string filename, Eigen::VectorXd u) {
    // Open file for data write
    std::fstream file;
    file.open(filename, std::ios::out);
    if (!file) {
        std::cerr << "Cant open file to save data: " << filename << std::endl;
        exit(1);
    }
    file << "u_n";
    file << std::endl;

    for (auto v : u) {
        file << v << "\n";
    }
    file.close();
}

auto main() -> int {
    // std::cout << "Hello N07" << std::endl;
    // solve_system_of_equations_newton();
    // solve_system_of_equations_damped();
    auto u = damped_newton_solve();
    save_to_file("./data/u_n.csv", u);
    return 0;
}