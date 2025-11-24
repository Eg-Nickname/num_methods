#include <chrono>
#include <fstream>
#include <iostream>
#include <ostream>

#include "bench.hpp"

long mesure_time_ms(std::function<void()> fn) {
  auto start = std::chrono::steady_clock::now();
  fn();
  auto end = std::chrono::steady_clock::now();
  return std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
      .count();
}

long mesure_time_micros(std::function<void()> fn) {
  auto start = std::chrono::steady_clock::now();
  fn();
  auto end = std::chrono::steady_clock::now();
  return std::chrono::duration_cast<std::chrono::microseconds>(end - start)
      .count();
}

long mesure_avg_time_ms(std::function<void()> fn, std::size_t runs) {
  long run_sum = 0;
  for (std::size_t i = 0; i < runs; i++) {
    run_sum += mesure_time_ms(fn);
  }
  return run_sum / runs;
}

long mesure_avg_time_micros(std::function<void()> fn, std::size_t runs) {
  long run_sum = 0;
  for (std::size_t i = 0; i < runs; i++) {
    run_sum += mesure_time_micros(fn);
  }
  return run_sum / runs;
}

void bench_method(BenchData bd) {
  // Open file for data write
  std::fstream file;
  file.open(bd.output_filename, std::ios::out);
  if (!file) {
    std::cerr << "Cant open file to save benchmark data: " << bd.output_filename
              << std::endl;
    exit(1);
  }
  file << bd.file_header;
  file << std::endl;

  for (long i = 0; i < bd.samples; i++) {
    long N =
        (bd.samples_upper_bound - bd.samples_lower_bound) / bd.samples * i +
        bd.samples_lower_bound;
    long time = mesure_time_micros([&bd, N]() { bd.fn(N); });
    if (bd.mesure_avg && time < bd.avg_ms_treshold * 1000) {
      time = mesure_avg_time_micros([&bd, N]() { bd.fn(N); }, 10);
    }
    file << N << ";" << time << ";"
         << std::endl; // endl to prevent data loss on crash or force strop
    std::cout << "For N=" << N << " sol_time = " << time << "\n";
  }
  file << std::endl;
  file.close();
}