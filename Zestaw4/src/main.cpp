#include "Eigen/Core"
#include <Eigen/Dense>
#include <iostream>

Eigen::Vector2d calc_eig_vals(Eigen::Matrix2d mat);

auto main() -> int {
  //   std::cout << "Hello nums! \n" << std::endl;
  Eigen::Matrix2d mat{{2, 2}, {2, 5}};
  //   Eigen::Vector2d ev = calc_eig_vals(mat);
  std::cout << "Eigen val of A mat is approx: "
            << calc_eig_vals(mat).transpose() << std::endl;
  Eigen::Matrix2d mat_a = mat - 6 * Eigen::Matrix2d::Identity();
  //   std::cout << mat_a << std::endl;

  std::cout << "Eigen val of mat a) is approx: "
            << calc_eig_vals(mat_a).transpose() << std::endl;

  std::cout << "Eigen val of mat b) is approx: "
            << calc_eig_vals(mat.inverse()).transpose() << std::endl;

  Eigen::Matrix2d mat_c = mat - 5 * Eigen::Matrix2d::Identity();
  std::cout << "Eigen val of mat c) is approx: "
            << calc_eig_vals(mat_c.inverse()).transpose() << std::endl;

  Eigen::Matrix2d mat_d = mat - 5.9999 * Eigen::Matrix2d::Identity();
  std::cout << "Eigen val of mat d) is approx: "
            << calc_eig_vals(mat_d.inverse()).transpose() << std::endl;

  Eigen::Matrix2d mat_e = mat - 0.9999 * Eigen::Matrix2d::Identity();
  std::cout << "Eigen val of mat e) is approx: "
            << calc_eig_vals(mat_e.inverse()).transpose() << std::endl;

  Eigen::Matrix2d mat_f = mat - 3.5 * Eigen::Matrix2d::Identity();
  std::cout << "Eigen val of mat f) is approx: "
            << calc_eig_vals(mat_f.inverse()).transpose() << std::endl;
  return 0;
}

Eigen::Vector2d calc_eig_vals(Eigen::Matrix2d mat) {
  Eigen::Vector2d ev;

  double eps = 1e-8;
  int max_iter = 1000;

  Eigen::Vector2d z_n = Eigen::Vector2d::Random();
  z_n.normalize();

  Eigen::Vector2d z_n1;
  int iter = 0;

  do {
    z_n1 = mat * z_n;
    z_n1.normalize();

    if (std::abs(z_n.dot(z_n1)) > 1.0 - eps)
      break;

    z_n = z_n1;
    iter++;
  } while (iter < max_iter);

  Eigen::Vector2d v1 = z_n1;
  ev(0) = v1.dot(mat * v1);

  Eigen::Vector2d z = Eigen::Vector2d::Random();
  z -= v1.dot(z) * v1;
  z.normalize();

  iter = 0;
  do {
    z_n1 = mat * z;

    // orthogonalize
    z_n1 -= v1.dot(z_n1) * v1;
    z_n1.normalize();

    if (std::abs(z.dot(z_n1)) > 1.0 - eps)
      break;

    z = z_n1;
    iter++;
  } while (iter < max_iter);

  Eigen::Vector2d v2 = z_n1;
  ev(1) = v2.dot(mat * v2);

  return ev;
}
