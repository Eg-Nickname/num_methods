#ifndef BENCH_HPP
#define BENCH_HPP
#include <functional>

long mesure_time_ms(std::function<void()> fn);
long mesure_time_micros(std::function<void()> fn);
long mesure_avg_time_ms(std::function<void()> fn, std::size_t runs);
long mesure_avg_time_micros(std::function<void()> fn, std::size_t runs);

#endif