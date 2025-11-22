#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <iostream>

#include "bench.hpp"
#include "solvers/dense_solver.hpp"
#include "solvers/shermanmorrison_solver.hpp"
#include "solvers/sparse_solver.hpp"
#include "vector_gen.hpp"

auto main() -> int {
  const long N = 1000;

  std::cout << "Sparse LU mat solve micros = "
            << mesure_time_micros([]() { solve_sp_mat_lu(N); }) << "\n";

  std::cout << "Sparse QR mat solve micros = "
            << mesure_time_micros([]() { solve_sp_mat_qr(N); }) << "\n";

  std::cout << "Sherman Morrison mat solve micros = "
            << mesure_avg_time_micros(
                   []() { solve_mat_sherman_morrison(N).head(5); }, 10)
            << "\n";

  return 0;
}
