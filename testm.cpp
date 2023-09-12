#include "matrix.hpp"
#include <iostream>
#include <random>

int main() {
  auto A = generateMatrix(10, 15, 100.0);
  std::stringstream SS;
  SS << A;
  auto B = Matrix<double>{SS.str()};
  if (A != B) {
    std::cerr << "Matrix should be same, but they arn't:" << std::endl;
    std::cerr << "A: \n" << A;
    std::cerr << "B: \n" << B;
    return -1;
  }
  std::cout << "OK" << std::endl;
  return 0;
}
