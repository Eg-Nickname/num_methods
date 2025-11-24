#ifndef BENCH_HPP
#define BENCH_HPP
#include <Eigen/Dense>
#include <functional>

struct BenchData {
  std::string output_filename;
  std::string file_header;
  std::function<void(long)> fn;
  long samples;
  long samples_lower_bound;
  long samples_upper_bound;
  bool mesure_avg;
  long avg_ms_treshold;
};

long mesure_time_ms(std::function<void()> fn);
long mesure_time_micros(std::function<void()> fn);
long mesure_avg_time_ms(std::function<void()> fn, std::size_t runs);
long mesure_avg_time_micros(std::function<void()> fn, std::size_t runs);

void bench_method(BenchData);

#endif