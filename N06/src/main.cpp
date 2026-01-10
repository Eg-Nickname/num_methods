#include "Eigen/Core"
#include <cmath>
#include <cstdlib>
#include <iostream>

using complex = std::complex<double>;

Eigen::VectorXcd calc_derivative_coef(const Eigen::VectorXcd &coef) {
    std::uint32_t coef_count = coef.size() - 1;
    Eigen::VectorXcd derivative_coef(coef_count);
    for (std::uint32_t i = 0; i < coef_count; i++) {
        derivative_coef(i) = (complex)(coef_count - i) * coef(i);
    }
    return derivative_coef;
}

Eigen::VectorXcd deflate(complex z, const Eigen::VectorXcd &coef) {
    int n = coef.size() - 1;
    Eigen::VectorXcd b(n);

    b(0) = coef(0);

    for (int i = 1; i < n; i++) {
        b(i) = coef(i) + z * b(i - 1);
    }

    return b;
}

complex horner_evaluate(complex z, const Eigen::VectorXcd &coef) {
    complex P = coef(0);
    for (int i = 1; i < coef.size(); i++) {
        P = (P * z) + coef(i);
    }
    return P;
}

class Laguerre {
    Eigen::VectorXcd coef;
    Eigen::VectorXcd first_deriv_coef;
    Eigen::VectorXcd second_deriv_coef;

    public:
    Laguerre(Eigen::VectorXcd coef) {
        this->coef = coef;
        this->first_deriv_coef = calc_derivative_coef(coef);
        this->second_deriv_coef = calc_derivative_coef(this->first_deriv_coef);
    }

    double get_n() { return (double)this->coef.size() - 1; }

    complex solve() { return solve(complex(0.1, 0.1)); }
    complex solve(complex guess) {
        complex z_n = guess;
        complex z_n1 = guess;
        double n = (double)this->coef.size() - 1;
        do {
            z_n = z_n1;
            complex p_z = horner_evaluate(z_n, this->coef);
            complex pd_z = horner_evaluate(z_n, this->first_deriv_coef);
            complex pdd_z = horner_evaluate(z_n, this->second_deriv_coef);

            complex sqrt = std::sqrt(
                (n - 1.0) * ((n - 1.0) * pd_z * pd_z - n * p_z * pdd_z));
            complex m1 = pd_z + sqrt;
            complex m2 = pd_z - sqrt;

            double m1_mag = std::abs(m1);
            double m2_mag = std::abs(m2);

            // Chose bigger denominator
            if (m1_mag < m2_mag) {
                z_n1 = z_n - ((n * p_z) / m2);
            } else {
                z_n1 = z_n - ((n * p_z) / m1);
            }
        } while (std::abs(z_n - z_n1) > 1e-8);
        return z_n1;
    }
};

// complex laguerre(const Eigen::VectorXcd &coef, complex starting) {}
std::vector<complex> lag(const Eigen::VectorXcd &coef) {
    std::vector<complex> sol{};
    Laguerre P = Laguerre(coef);
    sol.push_back(P.solve());
    // We dont need to smooth first
    Eigen::VectorXcd cur_coef = deflate(sol.back(), coef);

    while (cur_coef.size() > 3) {
        // Calculate root
        Laguerre P_n1 = Laguerre(cur_coef);
        auto z = P_n1.solve();
        // smooth it
        auto z_smooth = P.solve(z);
        sol.push_back(z_smooth);
        // deflate
        cur_coef = deflate(z_smooth, cur_coef);
    }
    // Calculate 2 remaning polynomial roots
    complex discriminant =
        std::sqrt(cur_coef(1) * cur_coef(1) - 4.0 * cur_coef(0) * cur_coef(2));

    complex z1 = (-cur_coef(1) + discriminant) / (2.0 * cur_coef(0));
    complex z2 = (-cur_coef(1) - discriminant) / (2.0 * cur_coef(0));

    // Smoothing of calculated roots
    sol.push_back(P.solve(z1));
    sol.push_back(P.solve(z2));

    return sol;
}
// (complex)(0)); }

auto main() -> int {
    Eigen::VectorXcd v = Eigen::VectorXcd::Zero(8);
    v(0) = 243;
    v(1) = -486;
    v(2) = 783;
    v(3) = -990;
    v(4) = 558;
    v(5) = -28;
    v(6) = -72;
    v(7) = 16;
    auto l = Laguerre(v);
    // std::cout << l.horner_evaluate((complex)(-2.4), v);
    // std::cout << l.solve();
    // auto deflated = deflate(l.solve(), v);
    auto sol = lag(v);
    for (auto c : sol) {
        std::cout << c << ", ";
    }
    std::cout << std::endl;
    // std::cout << horner_evaluate((complex)(2.0 / 3.0), deflated);

    return 0;
}
