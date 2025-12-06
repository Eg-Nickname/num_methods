#include "Eigen/Core"
#include <cassert>
#include <iostream>

#include "../typedefs.hpp"
#include "../vector_gen.hpp"
#include "thomas_solver.hpp"

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
  ThomasSolver(Eigen::VectorXd &&x, Eigen::VectorXd &&y, Eigen::VectorXd &&z)
      : y(std::move(y)), z(std::move(z)) {

    long N = this->y.size();
    this->c = Eigen::VectorXd(N);

    this->c(0) = x(0) / this->y(0);
    for (long i = 1; i < this->y.size(); i++) {
      this->c(i) = x(i) / (this->y(i) - this->z(i) * this->c(i - 1));
    }
  }

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

ThomasSolver create_solver_for_A(long N) {
  Eigen::VectorXd a = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd c = Eigen::VectorXd::Zero(N);

  float64_t h = 2 / ((float64_t)N - 1.0);

  for (long n = 0; n < N; n++) {
    a(n) = 1.0 / (h * h);
    b(n) = -2.0 / (h * h);
    c(n) = 1.0 / (h * h);
  }
  a(0) = 0;
  c(N - 1) = 0;
  b(0) = 1;
  b(N - 1) = 1;

  return ThomasSolver(std::move(a), std::move(b), std::move(c));
}

Eigen::VectorXd solve_mat_thomas(long N) {
  auto A_solver = create_solver_for_A(N);

  auto b = gen_b_vector(N);
  auto u = A_solver.solve(b);
  return u;
}