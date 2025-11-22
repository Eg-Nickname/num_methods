#include <cmath>

#include "typedefs.hpp"
#include "vector_gen.hpp"

Eigen::VectorXd gen_b_vector(long N) {
  Eigen::VectorXd b = Eigen::VectorXd(N);
  b.setZero();
  for (long n = 1; n <= N; n++) {
    float64_t val = std::cos((4.0 * std::numbers::pi * (n - 1)) / (double)N);
    b(n - 1) = val;
  }

  return b;
}