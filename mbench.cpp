#include "matrix.hpp"
#include <chrono>
#include <functional>
#include <iostream>
#include <map>

auto getCompilerId() {
#ifndef COMPILER_ID_NAME
#error "COMPILER_ID_NAME not set"
#endif
  return COMPILER_ID_NAME;
}

template <std::invocable Fn> auto bench(Fn &&F) {
  auto begin = std::chrono::high_resolution_clock::now();
  auto ret = std::invoke(std::forward<Fn>(F));
  auto end = std::chrono::high_resolution_clock::now();
  return std::make_tuple(std::move(ret), end - begin);
}

auto parseArgs(int argc, char **argv) {
  if (argc != 6) {
    throw std::runtime_error(
        "Wrong number of arguments. Expected: <dim from(log)> <dim to(log)> "
        "<thread from> <thread to> <extra points>");
  }
  auto dimFrom = std::stoul(argv[1]);
  auto dimTo = std::stoul(argv[2]);
  if (dimFrom > dimTo)
    throw std::runtime_error{"<dim to> must be not less <dim from>"};
  auto threadFrom = std::stoul(argv[3]);
  auto threadTo = std::stoul(argv[4]);
  auto extraPoints = std::stoul(argv[5]);
  if (extraPoints == 0)
    throw std::runtime_error{"<extra points> should be at least one"};
  if (threadFrom > threadTo)
    throw std::runtime_error{"<thread to> must be not less <thread from>"};
  return std::make_tuple(dimFrom, dimTo, threadFrom, threadTo, extraPoints);
}

auto runBench(unsigned dim, unsigned threads) {
  omp_set_num_threads(threads);
  auto a = generateMatrix(dim, dim, 10.0);
  auto b = generateMatrix(dim, dim, 10.0);

  auto &&[c, time] = bench([&]() { return Matrix<double>::multiply_st(a, b); });
  auto &&[c_par, time_par] =
      bench([&]() { return Matrix<double>::multiply_par(a, b); });
  auto &&[c_trans, time_trans] =
      bench([&]() { return Matrix<double>::multiply_st_transpose(a, b); });
  auto &&[c_par_trans, time_par_trans] =
      bench([&]() { return Matrix<double>::multiply_par_transpose(a, b); });

  if (c != c_par) {
    throw std::runtime_error{"par presented incorrect result"};
  }

  if (c != c_trans) {
    throw std::runtime_error{"st_trans presented incorrect result"};
  }

  if (c != c_par_trans) {
    throw std::runtime_error{"par_trans presented incorrect result"};
  }

  return std::make_tuple(time, time_par, time_trans, time_par_trans);
}

int main(int argc, char **argv) try {
  auto &&[dimFrom, dimTo, threadFrom, threadTo, extraPoints] =
      parseArgs(argc, argv);

  std::map<std::pair<unsigned, unsigned>,
           std::invoke_result_t<decltype(runBench), unsigned, unsigned>>
      results;
  for (auto dimLog = dimFrom; dimLog < dimTo + 1; ++dimLog) {
    for (auto extra = 0u; extra < extraPoints; ++extra) {
      auto dim = (1u << dimLog) + (1 << dimLog) * extra / extraPoints;
      for (auto nt = threadFrom; nt < threadTo + 1; ++nt)
        results[std::make_pair(dim, nt)] = runBench(dim, nt);
    }
  }
  for (auto &&[key, result] : results) {
    std::cout << getCompilerId() << ";" << key.first << "; " << key.second
              << "; " << std::get<0>(result).count() << ";"
              << std::get<1>(result).count() << ";"
              << std::get<2>(result).count() << ";"
              << std::get<3>(result).count() << ";" << std::endl;
  }
} catch (std::runtime_error &e) {
  std::cerr << e.what() << std::endl;
  return -1;
}
