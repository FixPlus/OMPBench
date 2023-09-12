#include "matrix.hpp"
#include <iostream>

auto parseArgs(int argc, char **argv) {
  if (argc != 4) {
    throw std::runtime_error(
        "Wrong number of arguments. Expected: <x> <y> <range>");
  }
  unsigned long x, y;
  x = std::stoul(argv[1]);
  y = std::stoul(argv[2]);
  double range;
  range = std::stod(argv[3]);
  if (range <= 0.0) {
    throw std::runtime_error{"range should be greater than zero"};
  }
  return std::make_tuple(x, y, range);
}

int main(int argc, char **argv) try {
  auto &&[x, y, range] = parseArgs(argc, argv);

  auto m = generateMatrix(x, y, range);
  std::cout << m;
  return 0;
} catch (std::runtime_error &e) {
  std::cerr << e.what() << std::endl;
  return -1;
}
