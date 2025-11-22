#ifndef DENSE_SOLVER_HPP
#define DENSE_SOLVER_HPP

#include <Eigen/Dense>

Eigen::MatrixXd gen_dense_A(long N);
Eigen::MatrixXd gen_dense_B(long N);

Eigen::VectorXd solve_d_mat_lu(long N);
Eigen::VectorXd solve_d_mat_qr(long N);
// Eigen::VectorXd solve_d_mat_LDLT(long N);

#endif