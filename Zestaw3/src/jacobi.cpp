#include <Eigen/Dense>
#include <iostream>

void Jacobi();
void GS();

int main() {

  // Jacobi();
  GS();
}

void Jacobi() {
  Eigen::MatrixXd A{{4, -1, 0}, {-1, 4, -1}, {0, -1, 4}};
  Eigen::Vector3d b{2, 6, 2};

  // Calculate inverse of diagonal matirx D and store it in column vector to
  // save space
  Eigen::Vector3d D_inv{};
  for (std::uint32_t i = 0; i < A.cols(); i++) {
    D_inv(i) = 1 / A(i, i);
  }
  // R matrix
  Eigen::MatrixXd R = A;
  R.diagonal().setZero();

  // Needed to calculate error rate in each itration
  Eigen::Vector3d x_analytic{1, 2, 1};

  Eigen::Vector3d x_n{0, 0, 0};
  for (std::size_t i = 0; i < 40; i++) {
    Eigen::VectorXd x_n1 = D_inv.asDiagonal().toDenseMatrix() * (b - R * x_n);
    std::cout << "Error at iter[" << i
              << "] = " << (x_analytic - x_n1).transpose() << std::endl;
    x_n = x_n1;
  }

  std::cout << "x_n after 50 iters = \n" << x_n << std::endl;
}

void GS() {
  Eigen::MatrixXd A{{4, -1, 0}, {-1, 4, -1}, {0, -1, 4}};
  Eigen::Vector3d b{2, 6, 2};

  Eigen::MatrixXd M(A.cols(), A.rows());
  Eigen::MatrixXd N(A.cols(), A.rows());

  for (std::uint32_t i = 0; i < A.rows(); i++) {
    for (std::uint32_t j = 0; j < A.cols(); j++) {
      if (i >= j) {
        M(i, j) = A(i, j);
        N(i, j) = 0;
      } else {
        M(i, j) = 0;
        N(i, j) = A(i, j);
      }
    }
  }
  std::cout << M << std::endl;
  std::cout << "N\n";
  std::cout << N << std::endl;
}