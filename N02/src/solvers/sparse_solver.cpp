#include <Eigen/Dense>

#include "../typedefs.hpp"
#include "../vector_gen.hpp"
#include "sparse_solver.hpp"

Eigen::SparseMatrix<double> gen_sparse_A(long N) {
  // matrix must be greater than 1
  float64_t h = 2 / ((float64_t)N - 1.0);

  Eigen::SparseMatrix<double> sp_mat = Eigen::SparseMatrix<double>(N, N);
  std::vector<Eigen::Triplet<double>> tripletList =
      std::vector<Eigen::Triplet<double>>(N);

  // for (int i = 0; i < N; i++) {
  //   tripletList.push_back(Eigen::Triplet<double>(i, i, -2.0 / (h * h)));
  //   if (i < N - 1)
  //     tripletList.push_back(Eigen::Triplet<double>(i, i + 1, 1.0 / (h * h)));
  //   if (i > 0)
  //     tripletList.push_back(Eigen::Triplet<double>(i, i - 1, 1.0 / (h * h)));
  // }

  for (long x = 1; x < N - 1; x++) {
    tripletList.push_back(Eigen::Triplet<double>(x, x, -2.0 / (h * h)));
    tripletList.push_back(Eigen::Triplet<double>(x, x + 1, 1.0 / (h * h)));
    tripletList.push_back(Eigen::Triplet<double>(x, x - 1, 1.0 / (h * h)));
  }

  // Diagonal fix
  tripletList.push_back(Eigen::Triplet<double>(0, 1, 1.0 / (h * h)));
  tripletList.push_back(Eigen::Triplet<double>(N - 1, N - 2, 1.0 / (h * h)));
  tripletList.push_back(Eigen::Triplet<double>(0, 0, -2.0 / (h * h)));
  tripletList.push_back(Eigen::Triplet<double>(N - 1, N - 1, -2.0 / (h * h)));

  // Corners
  tripletList.push_back(Eigen::Triplet<double>(N - 1, 0, 1.0 / (h * h)));
  tripletList.push_back(Eigen::Triplet<double>(0, N - 1, 1.0 / (h * h)));

  sp_mat.setFromTriplets(tripletList.begin(), tripletList.end());
  sp_mat.makeCompressed();

  return sp_mat;
}

// Based on
// https://libeigen.gitlab.io/eigen/docs-nightly/classEigen_1_1SparseLU.html
Eigen::VectorXd solve_sp_mat_lu(long N) {
  auto sp_mat = gen_sparse_A(N);
  auto b = gen_b_vector(N);

  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      sparse_lu_solver;
  sparse_lu_solver.analyzePattern(sp_mat);
  sparse_lu_solver.factorize(sp_mat);

  auto u = sparse_lu_solver.solve(b);
  return u;
}

Eigen::VectorXd solve_sp_mat_qr(long N) {
  auto sp_mat = gen_sparse_A(N);
  // Requirement
  // https://libeigen.gitlab.io/eigen/docs-nightly/classEigen_1_1SparseQR.html#a5f13e65437ada7e10a85de8cd55db11d
  sp_mat.makeCompressed();
  auto b = gen_b_vector(N);

  Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      sparse_qr_solver;
  sparse_qr_solver.analyzePattern(sp_mat);
  sparse_qr_solver.factorize(sp_mat);

  auto u = sparse_qr_solver.solve(b);
  return u;
}