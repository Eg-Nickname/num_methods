#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <cmath>
#include <iostream>

using float64_t = double;

bool thomas_solve(const Eigen::VectorXd &lower_in,
                  const Eigen::VectorXd &diag_in,
                  const Eigen::VectorXd &upper_in, Eigen::VectorXd &rhs,
                  double eps = 1e-14) {
  int n = diag_in.size();
  if ((int)lower_in.size() != n - 1 || (int)upper_in.size() != n - 1 ||
      (int)rhs.size() != n)
    return false;

  // make local copies of diagonals because we'll overwrite them
  Eigen::VectorXd c = upper_in; // modified super-diagonal
  Eigen::VectorXd d = diag_in;  // modified diagonal

  // Forward elimination
  for (int i = 1; i < n; ++i) {
    double pivot = d[i - 1];
    if (std::abs(pivot) < eps)
      return false; // pivot too small
    double m = lower_in[i - 1] / pivot;
    d[i] -= m * c[i - 1];
    rhs[i] -= m * rhs[i - 1];
  }

  // Back substitution
  double pivot = d[n - 1];
  if (std::abs(pivot) < eps)
    return false;
  rhs[n - 1] /= pivot;
  for (int i = n - 2; i >= 0; --i) {
    rhs[i] = (rhs[i] - c[i] * rhs[i + 1]) / d[i];
  }
  return true;
}

Eigen::VectorXd gen_N01_b_vector(long N) {
  Eigen::VectorXd b = Eigen::VectorXd(N);
  b.setZero();
  b(0) = 1;
  b(N - 1) = 1;
  return b;
}

Eigen::SparseMatrix<double> gen_sparse_A(long N) {
  // matrix must be greater than 1
  assert(N != 0);
  float64_t h = 2 / ((float64_t)N - 1.0);

  Eigen::SparseMatrix<double> sp_mat = Eigen::SparseMatrix<double>(N, N);
  std::vector<Eigen::Triplet<double>> tripletList =
      std::vector<Eigen::Triplet<double>>(N);

  for (int i = 0; i < N; i++) {
    tripletList.push_back(Eigen::Triplet<double>(i, i, -2.0 / (h * h)));
    if (i < N - 1)
      tripletList.push_back(Eigen::Triplet<double>(i, i + 1, 1.0 / (h * h)));
    if (i > 0)
      tripletList.push_back(Eigen::Triplet<double>(i, i - 1, 1.0 / (h * h)));
  }
  sp_mat.setFromTriplets(tripletList.begin(), tripletList.end());

  sp_mat.coeffRef(0, 0) = 1;
  sp_mat.coeffRef(N - 1, N - 1) = 1;

  return sp_mat;
}

Eigen::VectorXd solve_sp_mat_lu(long N) {
  auto sp_mat = gen_sparse_A(N);
  auto b = gen_N01_b_vector(N);

  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
      sparse_solver;
  sparse_solver.analyzePattern(sp_mat);
  sparse_solver.factorize(sp_mat);

  auto u = sparse_solver.solve(b);
  return u;
}

int main() {
  int N = 100000; // size
  double h = 2.0 / ((double)N - 1.0);
  Eigen::VectorXd lower = Eigen::VectorXd::Constant(N - 1, 1.0 / (h * h));
  Eigen::VectorXd diag = Eigen::VectorXd::Constant(N, -2.0 / (h * h));
  diag(0) = 1;
  diag(N - 1) = 1;
  Eigen::VectorXd upper = Eigen::VectorXd::Constant(N - 1, 1.0 / (h * h));

  Eigen::VectorXd b = gen_N01_b_vector(N);

  auto start = std::chrono::steady_clock::now();
  bool ok = thomas_solve(lower, diag, upper, b);
  auto end = std::chrono::steady_clock::now();
  if (!ok) {
    std::cerr << "Thomas failed: small pivot. Use a pivoting solver.\n";
    return 1;
  }
  auto t_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
                  .count();
  std::cout << "Thomas mat solve ms = " << t_ms << "\n";

  std::cout << "Thomas solution x (first 5 entries):\n" << b.head(5) << "\n";
  //   std::cout << b << std::endl;
  auto sp_start = std::chrono::steady_clock::now();
  auto u = solve_sp_mat_lu(N);
  auto sp_end = std::chrono::steady_clock::now();
  auto sp_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(sp_end - sp_start)
          .count();
  std::cout << "Sparse mat solve ms = " << sp_ms << "\n";

  std::cout << "Sparse LU solution u (first 5 entries):\n" << u.head(5) << "\n";

  return 0;
}
