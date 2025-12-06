#ifndef MATH_MATRIX_H
#define MATH_MATRIX_H

#pragma once

#include "Vector.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// Matrix
//
// an N⨉M matrix
//
template<int N, int M, typename T = double>
class Matrix
{
 public:
  using scalar_type = T;
  using eigen_type = Eigen::Matrix<T, N, M>;

 protected:
  eigen_type m_;

 public:
  // Construct an uninitialized Matrix.
  Matrix() = default;

  // Construct the Matrix `k I`.
  Matrix(T k = T{1}) requires (N == M)
  {
    m_.setIdentity();
    if (k != T{1})
      m_ *= k;
  }

  // Construct a Matrix from it eigen_type.
  Matrix(eigen_type const& m) : m_(m) { }

  template<typename Derived>
  Matrix(Eigen::MatrixBase<Derived> const& other) : m_(other) { }

  template<typename Derived>
  Matrix& operator=(Eigen::MatrixBase<Derived> const& other)
  {
    m_ = other;
    return *this;
  }

  eigen_type& eigen() { return m_; }
  eigen_type const& eigen() const { return m_; }

  T& operator()(int row, int col) { return m_(row, col); }
  T operator()(int row, int col) const { return m_(row, col); }

  // Return the determinant.
  T determinant() const requires (N == M) { return m_.determinant(); }

  bool isnan() const
  {
    bool has_nan = false;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        has_nan |= std::isnan(m_(i, j));
    return has_nan;
  }

  // Check if all coordinates are finite (not infinite or NaN).
  bool isfinite() const
  {
    bool finite = true;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        finite &= std::isfinite(m_(i, j));
    return finite;
  }

 public:
  // Add another matrix.
  Matrix& operator+=(Matrix const& m2)
  {
    m_ += m2.m_;
    return *this;
  }

  // Subtract another matrix.
  Matrix& operator-=(Matrix const& m2)
  {
    m_ -= m2.m_;
    return *this;
  }

  // Multiply the matrix with a scalar.
  Matrix& operator*=(T scalar)
  {
    m_ *= scalar;
    return *this;
  }

  // Divide the matrix by a scalar.
  Matrix& operator/=(T scalar)
  {
    m_ /= scalar;
    return *this;
  }

  // Multiply with another matrix.
  template<int P>
  Matrix<N, P, T> operator*(Matrix<M, P, T> const& m2) const
  {
    return {m_ * m2.eigen()};
  }

  // Return the inverse of the matrix.
  Matrix inverse() const requires (N == M) { return {m_.inverse()}; }

  // Return the transpose of the matrix.
  Matrix<M, N, T> transpose() const { return {m_.transpose()}; }

  // Negate this matrix.
  void negate()
  {
    m_ = -m_;
  }

  // Return the negated matrix.
  Matrix operator-() const
  {
    return {-m_};
  }

  // Add another matrix.
  Matrix operator+(Matrix const& m2) const
  {
    return {m_ + m2.m_};
  }

  // Subtract another matrix.
  Matrix operator-(Matrix const& m2) const
  {
    return {m_ - m2.m_};
  }

  // Multiply by a scalar.
  Matrix operator*(T scalar) const
  {
    return {m_ * scalar};
  }

  // Divide by a scalar.
  Matrix operator/(T scalar) const
  {
    return {m_ / scalar};
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << '[';
    for (int i = 0; i < N; ++i)
    {
      char const* separator = "";
      for (int j = 0; j < M; ++j)
      {
        os << separator << m_(i, j);
        separator = ", ";
      }
      if (i != N - 1)
        os << '\n';
    }
    os << ']';
  }
#endif
};

template<int N, int M, typename T>
Vector<N, T> operator*(Matrix<N, M, T> const& m, Vector<M, T> const& v)
{
  return {m.eigen() * v.eigen()};
}

template<int N, int M, typename T>
Vector<M, T> operator*(Vector<N, T> const& v, Matrix<N, M, T> const& m)
{
  return {v.eigen().transpose() * m.eigen()};
}

template<int N, int M, typename T>
Matrix<N, M, T> operator*(T scalar, Matrix<N, M, T> const& m)
{
  return {scalar * m.eigen()};
}

} // namespace math

#endif // MATH_MATRIX_H
