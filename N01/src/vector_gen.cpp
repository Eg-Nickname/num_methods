#include <cmath>

#include "typedefs.hpp"
#include "vector_gen.hpp"

Eigen::VectorXd gen_b_vector(long N) {
  Eigen::VectorXd b = Eigen::VectorXd(N);
  b.setZero();
  b(0) = 1;
  b(N - 1) = 1;
  return b;
}