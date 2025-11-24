#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cassert>
#include <exception>
#include <iostream>
#include <string>

#include "bench.hpp"
#include "solvers/dense_solver.hpp"
#include "solvers/shermanmorrison_solver.hpp"
#include "solvers/sparse_solver.hpp"
#include "vector_gen.hpp"

BenchData parse_b_flag(int argc, char **argv);
void parse_p_flag(int argc, char **argv);
long parse_num(std::string num_str, long err_val = 1);

auto main(int argc, char **argv) -> int {
  if (argc == 1) {
    std::cout << (gen_dense_A(200) * solve_d_mat_fullpiv_lu(200)).head(5)
              << "\n";
    std::cout << (gen_dense_A(200) * solve_d_mat_fullpiv_qr(200)).head(5)
              << "\n";

    std::cout << "exectute program with -b or -p flag for more info";
    return 0;
  } else if (std::string(argv[1]) == "-b" && argc < 7) {

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
  } else if (std::string(argv[1]) == "-p") {
    // Do usage of -p flag
  }

  if (std::string(argv[1]) == "-b") {
    parse_b_flag(argc, argv);
    bench_method(parse_b_flag(argc, argv));
  }

  if (std::string(argv[1]) == "-p") {
    // parse_p_flag(argc, argv);
  }

  return 0;
}

BenchData parse_b_flag(int argc, char **argv) {
  BenchData bd;
  bd.output_filename = argv[2];
  std::string method = argv[3];

  if (method == "SPARSE_LU") {
    bd.file_header = "n;Sparse LU;";
    bd.fn = [](long N) { solve_sp_mat_lu(N); };
  } else if (method == "SPARSE_QR") {
    bd.file_header = "n;Sparse QR;";
    bd.fn = [](long N) { solve_sp_mat_qr(N); };
  } else if (method == "DENSE_FULL_LU") {
    bd.file_header = "n;Dense Full Pivot LU;";
    bd.fn = [](long N) { solve_d_mat_fullpiv_lu(N); };
  } else if (method == "DENSE_PAR_LU") {
    bd.file_header = "n;Dense Partial Pivot LU;";
    bd.fn = [](long N) { solve_d_mat_partialpiv_lu(N); };
  } else if (method == "DENSE_FULL_QR") {
    bd.file_header = "n;Dense Full Pivot QR;";
    bd.fn = [](long N) { solve_d_mat_fullpiv_qr(N); };
  } else if (method == "n;DENSE_PAR_QR;") {
    bd.file_header = "n;Dense QR;";
    bd.fn = [](long N) { solve_d_mat_partialpiv_qr(N); };
  } else if (method == "SHERMAN_MORRISON") {
    bd.file_header = "n;Sherman Morrison;";
    bd.fn = [](long N) { solve_mat_sherman_morrison(N); };
  }

  bd.samples = parse_num(argv[4], 10);
  bd.samples_lower_bound = parse_num(argv[5], 10);
  bd.samples_upper_bound = parse_num(argv[6], 1000);
  bd.mesure_avg = false;
  bd.avg_ms_treshold = 10;
  if (argc == 9) {
    bd.mesure_avg = parse_num(argv[7], 0) == 1;
    bd.avg_ms_treshold = parse_num(argv[8], 20);
  }
  return bd;
}

long parse_num(std::string num_str, long err_val) {
  try {
    return std::stol(num_str);
  } catch (std::exception _) {
    std::cout << "Error parsing: " << num_str << "\n";
    return err_val;
  }
}
