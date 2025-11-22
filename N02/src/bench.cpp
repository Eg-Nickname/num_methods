#include <chrono>

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