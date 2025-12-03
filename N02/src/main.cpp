#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <exception>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>

#include "Eigen/Core"
#include "bench.hpp"
#include "solvers/dense_solver.hpp"
#include "solvers/shermanmorrison_solver.hpp"
#include "solvers/sparse_solver.hpp"
#include "typedefs.hpp"
#include "vector_gen.hpp"

BenchData parse_b_flag(int argc, char **argv);
void parse_p_flag(int argc, char **argv);
long parse_num(std::string num_str, long err_val = 1);

void gen_p_flag_data(char **argv);

auto main(int argc, char **argv) -> int {
  if (argc == 1) {
    std::cout << (gen_dense_A(200) * solve_d_mat_fullpiv_lu(200)).head(5)
              << "\n";
    std::cout << (gen_dense_A(200) * solve_d_mat_fullpiv_qr(200)).head(5)
              << "\n";

    std::cout << "exectute program with -b or -p flag for more info";
    return 0;
  } else if (std::string(argv[1]) == "-b" && argc < 8) {

    std::cout
        << "Execute program: " << argv[0]
        << "-b OUTPUT_FILE METHOD SAMPLES SAMPLE_RANGE_LOWER_BOUND "
           "SAMPLE_RANGE_UPPER_BOUND (OPTIONAL: MESURE_AVG AVG_MS_TRESHHOLD)"
        << std::endl;

    std::cout << "Methods:\n \t - SPARSE_LU\n \t - SPARSE_QR \n \t - "
                 "DENSE_FULL_LU \n \t - DENSE_FULL_QR\n "
                 "\t - DENSE_PAR_LU\n \t - DENSE_PAR_QR\n \t - SHERMAN_MORRISON"
              << std::endl;
    return 0;
  } else if (std::string(argv[1]) == "-p" && argc < 2) {
    // Do usage of -p flag
    return 0;
  }

  if (std::string(argv[1]) == "-b") {
    parse_b_flag(argc, argv);
    bench_method(parse_b_flag(argc, argv));
  }

  if (std::string(argv[1]) == "-p") {
    // parse_p_flag(argc, argv);
    gen_p_flag_data(argv);
  }

  return 0;
}

BenchData parse_b_flag(int argc, char **argv) {
  BenchData bd;
  bd.output_filename = argv[2];
  std::string method = argv[3];
  std::string col_name = argv[4];
  bd.file_header = std::string("n;") + col_name + ";";
  if (method == "SPARSE_LU") {
    bd.fn = [](long N) { solve_sp_mat_lu(N); };
  } else if (method == "SPARSE_QR") {
    bd.fn = [](long N) { solve_sp_mat_qr(N); };
  } else if (method == "DENSE_FULL_LU") {
    bd.fn = [](long N) { solve_d_mat_fullpiv_lu(N); };
  } else if (method == "DENSE_PAR_LU") {
    bd.fn = [](long N) { solve_d_mat_partialpiv_lu(N); };
  } else if (method == "DENSE_FULL_QR") {
    bd.fn = [](long N) { solve_d_mat_fullpiv_qr(N); };
  } else if (method == "DENSE_PAR_QR") {
    std::cout << "PARQR" << std::endl;
    bd.fn = [](long N) { solve_d_mat_partialpiv_qr(N); };
  } else if (method == "SHERMAN_MORRISON") {
    bd.fn = [](long N) { solve_mat_sherman_morrison(N); };
  } else {
    std::cout << "NO FN SELECTED" << std::endl;
  }

  bd.samples = parse_num(argv[5], 10);
  bd.samples_lower_bound = parse_num(argv[6], 10);
  bd.samples_upper_bound = parse_num(argv[7], 1000);
  bd.mesure_avg = false;
  bd.avg_ms_treshold = 10;
  if (argc == 10) {
    bd.mesure_avg = parse_num(argv[8], 0) == 1;
    bd.avg_ms_treshold = parse_num(argv[9], 20);
  }
  return bd;
}

long parse_num(std::string num_str, long err_val) {
  try {
    return std::stol(num_str);
  } catch (std::exception &_) {
    std::cout << "Error parsing: " << num_str << "\n";
    return err_val;
  }
}

void gen_p_flag_data(char **argv) {
  std::fstream file;
  file.open("./bench_data/solution_data.csv", std::ios::out);
  if (!file) {
    std::cerr << "Cant open file to save benchmark data: "
                 "./bench_data/solution_data.csv"
              << std::endl;
    exit(1);
  }
  file << "x_n;y_n";
  file << std::endl;
  const long N = parse_num(argv[2], 1000);
  float64_t h = 2.0 / (float64_t)(N - 1);
  Eigen::VectorXd u = solve_mat_sherman_morrison(N);
  for (int i = 0; i < N; i++) {
    // save data to file
    float64_t x = i * h;
    float64_t y = u(i);
    file << x << ";" << y << "\n";
  }
  file << std::endl;
}