#include <algorithm>
#include <cassert>
#include <concepts>
#include <fstream>
#include <random>
#include <ranges>
#include <span>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <omp.h>

template <typename T>
    requires std::integral<T> || std::floating_point<T> class Matrix {
public:
  friend class RowIterator;

  template <typename C> class RowIterator {
  public:
    using owner_type =
        std::conditional_t<std::is_const_v<C>, const Matrix, Matrix>;
    using value_type = std::span<C>;
    using reference = std::span<C>;
    using pointer = void;
    using difference_type = long;
    using iterator_category = std::input_iterator_tag;
    RowIterator(owner_type &matrix, unsigned long row)
        : m_matrix{&matrix}, m_row{row} {};

    auto operator*() const {
      return std::span<C>{m_matrix->m_data.data() + m_row * m_matrix->m_dimX,
                          m_matrix->m_dimX};
    }

    auto &operator+=(long Diff) {
      m_row += Diff;
      return *this;
    }

    auto operator[](unsigned long Id) const { return *(*this + Id); }

    auto &operator-=(long Diff) {
      m_row -= Diff;
      return *this;
    }

    auto operator+(long Diff) const {
      auto New = *this;
      New += Diff;
      return New;
    }

    auto operator-(long Diff) const { return *this + (-Diff); }

    auto operator<(const RowIterator &another) const {
      assert(m_matrix == another.m_matrix);
      return m_row < another.m_row;
    }
    auto operator>(const RowIterator &another) const {
      return another < *this && *this != another;
    }

    auto operator>=(const RowIterator &another) const {
      return !(*this < another);
    }

    auto operator<=(const RowIterator &another) const {
      return !(*this > another);
    }

    auto operator-(const RowIterator &another) const {
      assert(m_matrix == another.m_matrix);
      return m_row - another.m_row;
    }

    auto &operator++() {
      m_row++;
      return *this;
    }

    auto &operator++(int) {
      auto Old = *this;
      m_row++;
      return Old;
    }

    bool operator==(const RowIterator &another) const {
      return (m_matrix == another.m_matrix) && m_row == another.m_row;
    }
    bool operator!=(const RowIterator &another) const {
      return !(*this == another);
    }

  private:
    owner_type *m_matrix;
    unsigned long m_row;
  };
  Matrix(unsigned long x, unsigned long y)
      : m_dimX{x}, m_dimY{y}, m_data(
                                  [&]() {
                                    m_checkDimensions(x, y);
                                    return x * y;
                                  }(),
                                  0) {}
  Matrix(std::string Encoded) {
    std::istringstream IS{Encoded};
    IS >> m_dimX;
    if (IS.bad())
      throw std::runtime_error(
          "Failed to read matrix from string. Could not read dimX");
    IS >> m_dimY;
    if (IS.bad())
      throw std::runtime_error(
          "Failed to read matrix from string. Could not read dimY");
    m_checkDimensions(m_dimX, m_dimY);
    m_data.reserve(m_dimX * m_dimY);
    for (auto i = 0u; i < m_dimY; ++i)
      for (auto j = 0u; j < m_dimX; ++j) {
        T Val;
        IS >> Val;
        if (IS.bad())
          throw std::runtime_error(
              "Failed to read (" + std::to_string(j) + ", " +
              std::to_string(i) + ") value for matrix with dimensions of (" +
              std::to_string(m_dimX) + ", " + std::to_string(m_dimY) + ")");

        m_data.emplace_back(std::move(Val));
      }
    assert(m_data.size() == m_dimX * m_dimY);
  }

  auto dimX() const { return m_dimX; }
  auto dimY() const { return m_dimY; }

  auto &at(unsigned long x, unsigned long y) const {
    if (x >= m_dimX || y >= m_dimY) {
      std::ostringstream OS;
      OS << "(" << x << ", " << y << ") is out of range for matrix (" << m_dimX
         << ", " << m_dimY << ")";
      throw std::out_of_range{OS.str()};
    }
    return m_data.at(y * m_dimX + x);
  }

  auto &at(unsigned long x, unsigned long y) {
    if (x >= m_dimX || y >= m_dimY) {
      std::ostringstream OS;
      OS << "(" << x << ", " << y << ") is out of range for matrix (" << m_dimX
         << ", " << m_dimY << ")";
      throw std::out_of_range{OS.str()};
    }
    return m_data.at(y * m_dimX + x);
  }

  auto begin() { return RowIterator<T>{*this, 0u}; }

  auto begin() const { return RowIterator<const T>{*this, 0u}; }

  auto end() { return RowIterator<T>{*this, m_dimY}; }

  auto end() const { return RowIterator<const T>{*this, m_dimY}; }

  template <std::equality_comparable_with<T> P>
  bool operator==(const Matrix<P> &another) const {
    if (m_dimX != another.dimX() || m_dimY != another.dimY())
      return false;
    for (auto i = 0u; i < m_dimY; ++i)
      for (auto j = 0u; j < m_dimX; ++j) {
        auto &&a = at(j, i);
        auto &&b = another.at(j, i);
        if constexpr (std::floating_point<T> || std::floating_point<P>) {
          if (std::abs(a - b) >= FPrecision)
            return false;
        } else if (a != b)
          return false;
      }
    return true;
  }

  bool operator!=(const Matrix &another) const { return !(*this == another); }

  static Matrix loadFromFile(std::string_view filename) {
    std::ifstream file{filename};
    if (!file.is_open())
      throw std::runtime_error{"Cannot open file '" + std::string{filename} +
                               "' for reading"};
    file.seekg(std::ios::end);
    auto filesize = file.tellg();
    file.seekg(std::ios::beg);
    std::string buffer;
    buffer.resize(filesize);
    file.read(buffer.data(), filesize);
    return Matrix{buffer};
  }

  static void checkDimForMultiply(const Matrix &a, const Matrix &b) {
    if (a.dimX() == b.dimY())
      return;
    std::stringstream SS;
    SS << "Matrices a(" << a.dimX() << ", " << a.dimY() << ") and ";
    SS << "b(" << b.dimX() << ", " << b.dimY() << ")";
    SS << " cannot be multiplied: ax != by";
    throw std::runtime_error(SS.str());
  }

  static Matrix multiply_st(const Matrix &a, const Matrix &b) {
    checkDimForMultiply(a, b);
    Matrix Ret{b.dimX(), a.dimY()};
    auto &arrA = a.m_data;
    auto &arrB = b.m_data;
    auto &arrRet = Ret.m_data;
    const auto dimX = Ret.dimX();
    const auto dimY = Ret.dimY();
    const auto dimK = b.dimY();

    for (unsigned j = 0u; j < dimY; j++)
      for (unsigned i = 0u; i < dimX; i++) {
        arrRet[j * dimX + i] = 0u;
        for (int k = 0ul; k < dimK; ++k)
          arrRet[j * dimX + i] += arrA[k + j * dimK] * arrB[i + k * dimX];
      }

    return Ret;
  }

  static Matrix multiply_st_transpose(const Matrix &a, const Matrix &b) {
    checkDimForMultiply(a, b);
    Matrix Ret{b.dimX(), a.dimY()};
    auto bTrans = b.transpose();
    auto &arrA = a.m_data;
    auto &arrB = bTrans.m_data;
    auto &arrRet = Ret.m_data;
    const auto dimX = Ret.dimX();
    const auto dimY = Ret.dimY();
    const auto dimK = b.dimY();

    for (unsigned j = 0u; j < dimY; j++)
      for (unsigned i = 0u; i < dimX; i++) {
        arrRet[j * dimX + i] = 0u;
        for (int k = 0ul; k < dimK; ++k)
          arrRet[j * dimX + i] += arrA[k + j * dimK] * arrB[k + i * dimK];
      }

    return Ret;
  }
  static Matrix multiply_par(const Matrix &a, const Matrix &b) {
    checkDimForMultiply(a, b);
    Matrix Ret{b.dimX(), a.dimY()};
    auto &arrA = a.m_data;
    auto &arrB = b.m_data;
    auto &arrRet = Ret.m_data;
    const auto dimX = Ret.dimX();
    const auto dimY = Ret.dimY();
    const auto dimK = b.dimY();

#pragma omp parallel for collapse(2) firstprivate(dimX, dimY, dimK),           \
    shared(arrRet, arrA, arrB)
    for (unsigned j = 0u; j < dimY; j++)
      for (unsigned i = 0u; i < dimX; i++) {
        arrRet[j * dimX + i] = 0u;
        for (int k = 0ul; k < dimK; ++k)
          arrRet[j * dimX + i] += arrA[k + j * dimK] * arrB[i + k * dimX];
      }

    return Ret;
  }

  static Matrix multiply_par_transpose(const Matrix &a, const Matrix &b) {
    checkDimForMultiply(a, b);
    Matrix Ret{b.dimX(), a.dimY()};
    auto bTrans = b.transpose();
    auto &arrA = a.m_data;
    auto &arrB = bTrans.m_data;
    auto &arrRet = Ret.m_data;
    const auto dimX = Ret.dimX();
    const auto dimY = Ret.dimY();
    const auto dimK = b.dimY();

#pragma omp parallel for collapse(2) firstprivate(dimX, dimY, dimK),           \
    shared(arrRet, arrA, arrB)
    for (unsigned j = 0u; j < dimY; j++)
      for (unsigned i = 0u; i < dimX; i++) {
        arrRet[j * dimX + i] = 0u;
        for (int k = 0ul; k < dimK; ++k)
          arrRet[j * dimX + i] += arrA[k + j * dimK] * arrB[k + i * dimK];
      }

    return Ret;
  }
  Matrix transpose() const {
    Matrix Ret{dimY(), dimX()};
    auto &arrA = m_data;
    auto &arrRet = Ret.m_data;
    for (auto j = 0u; j < dimY(); ++j)
      for (auto i = 0u; i < dimX(); ++i)
        arrRet[i * dimY() + j] = arrA[j * dimX() + i];
    return Ret;
  }

  static constexpr auto FPrecision = 0.0001;
  static constexpr auto MaxDim = 1u << 15;

private:
  static void m_checkDimensions(unsigned long x, unsigned long y) {
    if (x > MaxDim)
      throw std::runtime_error(
          "Invalid matrix dimension x: " + std::to_string(x) +
          ". Should be less than " + std::to_string(MaxDim));
    if (x == 0)
      throw std::runtime_error("Invalid matrix dimension x: " +
                               std::to_string(x) + ". Should be at least 1");
    if (y == 0)
      throw std::runtime_error("Invalid matrix dimension y: " +
                               std::to_string(y) + ". Should be at least 1");
    if (y > MaxDim)
      throw std::runtime_error(
          "Invalid matrix dimension y: " + std::to_string(y) +
          ". Should be less than " + std::to_string(MaxDim));
  }
  unsigned long m_dimX;
  unsigned long m_dimY;
  std::vector<T> m_data;
};

template <typename T>
std::ostream &operator<<(std::ostream &OS, const Matrix<T> &matrix) {
  OS << matrix.dimX() << " " << matrix.dimY() << " " << std::endl;
  for (auto &&row : matrix) {
    for (auto &&Elem : row)
      OS << Elem << " ";
    OS << std::endl;
  }
  return OS;
}

auto &randomEngine() {
  static std::random_device R;
  static std::default_random_engine E(R());
  return E;
}

template <std::integral T> T genInRange(T From, T To) {
  assert(To > From && "At least one element in range expected");
  std::uniform_int_distribution<T> uniform_dist(From, To - 1);
  return uniform_dist(randomEngine());
}

template <typename T> T genInRange(T From, T To) {
  assert(To > From && "At least one element in range expected");
  std::uniform_real_distribution<T> uniform_dist(From, To - 1);
  return uniform_dist(randomEngine());
}

template <typename T>
auto generateMatrix(unsigned long x, unsigned y, T range) {
  Matrix<T> m{x, y};
  for (auto &&row : m)
    for (auto &e : row)
      e = genInRange(-range, range);
  return m;
}
