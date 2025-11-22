#ifndef SPARSE_SOLVER_HPP
#define SPARSE_SOLVER_HPP
#include <Eigen/Sparse>

Eigen::SparseMatrix<double> gen_sparse_A(long N);

Eigen::VectorXd solve_sp_mat_lu(long N);
Eigen::VectorXd solve_sp_mat_qr(long N);

#endif