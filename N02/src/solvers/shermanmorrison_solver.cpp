#include "Eigen/Core"
#include <cassert>

#include "../typedefs.hpp"
#include "../vector_gen.hpp"
#include "shermanmorrison_solver.hpp"

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

// Create specific solver for our B matrix (A = B + uuT)
ThomasSolver create_solver_for_B(long N, double tr, double bl) {
  Eigen::VectorXd a = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
  Eigen::VectorXd c = Eigen::VectorXd::Zero(N);

  float64_t h = 2 / ((float64_t)N - 1.0);

  for (long n = 0; n < N; n++) {
    a(n) = 1.0 / (h * h);
    b(n) = -2.0 / (h * h);
    c(n) = 1.0 / (h * h);
  }
  // Fix B_(1,1) and B_(N, N)
  float64_t gamma = -b(0);
  b(0) -= gamma;
  b(N - 1) -= (tr * bl) / gamma;

  return ThomasSolver(std::move(a), std::move(b), std::move(c));
}

Eigen::VectorXd solve_mat_sherman_morrison(long N) {
  // Generate B matrix as threediagonal
  // Use thomas algorithm to solve intermidied equations for sherman morrison
  // solve for u using Sherman Morrison

  // A = B + uvT
  float64_t h = 2.0 / (N - 1);
  float64_t gamma = 2.0 / (h * h);

  // top left corner val
  float64_t tr = 1.0 / (h * h);
  // bototm lef corner val
  float64_t bl = 1.0 / (h * h);

  Eigen::VectorXd u = Eigen::VectorXd::Zero(N);
  u(0) = gamma;
  u(N - 1) = 1.0 / (h * h);
  Eigen::VectorXd v = Eigen::VectorXd::Zero(N);
  v(0) = 1;
  v(N - 1) = (bl * tr) / gamma;

  auto B = create_solver_for_B(N, tr, bl);

  auto b = gen_b_vector(N);
  // y = B^-1 b => By = b
  auto y = B.solve(b);
  // z = B^-1 u = Bz = u
  auto z = B.solve(u);

  float64_t p = (v.transpose() * z);
  Eigen::VectorXd x = y - (z * (v.transpose() * y)) / (1.0 + p);

  return x;
}