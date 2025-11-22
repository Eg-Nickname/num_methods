#include <iostream>

#include "../typedefs.hpp"
#include "../vector_gen.hpp"
#include "dense_solver.hpp"

Eigen::MatrixXd gen_dense_B(long N) {
  // matrix must be greater than 1
  assert(N != 0);
  float64_t h = 2 / ((float64_t)N - 1.0);
  Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(N, N);

  for (long x = 1; x < N - 1; x++) {
    mat(x, x - 1) = 1.0 / (h * h);
    mat(x, x + 1) = 1.0 / (h * h);
    mat(x, x) = -2.0 / (h * h);
  }
  mat(0, 0) = -3.0 / (h * h);
  mat(N - 1, N - 1) = -3.0 / (h * h);

  // Diagonal fix
  mat(0, 1) = 1.0 / (h * h);
  mat(N - 1, N - 2) = 1.0 / (h * h);

  return mat;
}

// A = B + uuT
Eigen::MatrixXd gen_dense_A(long N) {
  auto mat = gen_dense_B(N);
  float64_t h = 2 / ((float64_t)N - 1.0);

  // Distortion from 3 diagonal
  mat(0, N - 1) = 1.0 / (h * h);
  mat(N - 1, 0) = 1.0 / (h * h);

  mat(0, 0) += 1.0 / (h * h);
  mat(N - 1, N - 1) += 1.0 / (h * h);
  return mat;
}

Eigen::VectorXd solve_d_mat_lu(long N) {
  auto d_mat = gen_dense_A(N);
  auto b = gen_b_vector(N);

  auto dense_lu = d_mat.fullPivLu();

  auto u = dense_lu.solve(b);
  return u;
}

Eigen::VectorXd solve_d_mat_qr(long N) {
  auto d_mat = gen_dense_A(N);
  auto b = gen_b_vector(N);

  auto dense_qr = d_mat.householderQr();

  auto u = dense_qr.solve(b);
  return u;
}
