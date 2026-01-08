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

template<int N, int M, typename T>
class MatrixData
{
 public:
  using eigen_type = Eigen::Matrix<T, N, M>;

 protected:
  eigen_type m_;

 public:
  // Construct an uninitialized Matrix.
  MatrixData() = default;

  // Construct the Matrix `k I`.
  explicit MatrixData(T k) requires (N == M)
  {
    m_.setIdentity();
    if (k != T{1})
      m_ *= k;
  }

  // Construct a Matrix from its eigen_type.
  MatrixData(eigen_type const& m) : m_(m) { }

  template<typename Derived>
  MatrixData(Eigen::MatrixBase<Derived> const& other) : m_(other) { }

  template<typename Derived>
  MatrixData& operator=(Eigen::MatrixBase<Derived> const& other)
  {
    m_ = other;
    return *this;
  }

  // Accessors.
  eigen_type& eigen() { return m_; }
  eigen_type const& eigen() const { return m_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

template<typename DerivedTypes>
struct MatrixOps
{
 private:
  static constexpr int derived_rows = DerivedTypes::rows;
  static constexpr int derived_cols = DerivedTypes::cols;
  using derived_scalar_type    = typename DerivedTypes::scalar_type;
  using derived_type           = typename DerivedTypes::derived_type;
  template<int R, int C>
  using derived_matrix_type    = typename DerivedTypes::template matrix_type<R, C>;
  template<int D>
  using derived_vector_type    = typename DerivedTypes::template vector_type<D>;

  auto& raw_() { return static_cast<derived_type*>(this)->raw(); }
  auto const& raw_() const { return static_cast<derived_type const*>(this)->raw(); }

 public:
  derived_scalar_type& operator()(int row, int col) { return raw_().operator()(row, col); }
  derived_scalar_type operator()(int row, int col) const { return raw_().operator()(row, col); }

  derived_scalar_type determinant() const requires (derived_rows == derived_cols) { return raw_().determinant(); }
  bool isnan() const { return raw_().isnan(); }
  bool isfinite() const { return raw_().isfinite(); }

  derived_type& operator+=(derived_type const& m2) { raw_().operator+=(m2.raw()); return static_cast<derived_type&>(*this); }
  derived_type& operator-=(derived_type const& m2) { raw_().operator-=(m2.raw()); return static_cast<derived_type&>(*this); }
  derived_type& operator*=(derived_scalar_type scalar) { raw_().operator*=(scalar); return static_cast<derived_type&>(*this); }
  derived_type& operator/=(derived_scalar_type scalar) { raw_().operator/=(scalar); return static_cast<derived_type&>(*this); }

  template<int P>
  derived_matrix_type<derived_rows, P> operator*(derived_matrix_type<derived_cols, P> const& m2) const
  {
    return static_cast<derived_matrix_type<derived_rows, P>>(raw_().operator*(m2.raw()));
  }

  derived_vector_type<derived_rows> operator*(derived_vector_type<derived_cols> const& v) const
  {
    return static_cast<derived_vector_type<derived_rows>>(raw_().operator*(v.raw()));
  }

  derived_type inverse() const requires (derived_rows == derived_cols) { return static_cast<derived_type>(raw_().inverse()); }
  derived_matrix_type<derived_cols, derived_rows> transpose() const { return static_cast<derived_matrix_type<derived_cols, derived_rows>>(raw_().transpose()); }

  void negate() { raw_().negate(); }
  derived_type operator-() const { return static_cast<derived_type>(raw_().operator-()); }

  derived_type operator+(derived_type const& m2) const { return static_cast<derived_type>(raw_().operator+(m2.raw())); }
  derived_type operator-(derived_type const& m2) const { return static_cast<derived_type>(raw_().operator-(m2.raw())); }
  derived_type operator*(derived_scalar_type scalar) const { return static_cast<derived_type>(raw_().operator*(scalar)); }
  derived_type operator/(derived_scalar_type scalar) const { return static_cast<derived_type>(raw_().operator/(scalar)); }
};

// Forward declaration, required for MatrixTypes<N, T>::derived_type.
template<int N, int M, typename T>
class Matrix;

template<int N, int M, typename T>
struct MatrixTypes
{
  static constexpr int rows = N;
  static constexpr int cols = M;
  using scalar_type    = T;
  using derived_type   = Matrix<N, M, T>;

  template<int R, int C>
  using matrix_type    = Matrix<R, C, T>;

  template<int D>
  using vector_type    = Vector<D, T>;
};

// Specialization of the Matrix operators specifically for math::Matrix itself.
template<int N, int M, typename T>
struct MatrixOps<MatrixTypes<N, M, T>>
{
 private:
  // Get underlying Eigen type.
  auto& eigen_() { return static_cast<Matrix<N, M, T>*>(this)->eigen(); }
  auto const& eigen_() const { return static_cast<Matrix<N, M, T> const*>(this)->eigen(); }

 public:
  T& operator()(int row, int col) { return eigen_()(row, col); }
  T operator()(int row, int col) const { return eigen_()(row, col); }

  T determinant() const requires (N == M) { return eigen_().determinant(); }

  bool isnan() const
  {
    bool has_nan = false;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        has_nan |= std::isnan(eigen_()(i, j));
    return has_nan;
  }

  // Check if all coordinates are finite (not infinite or NaN).
  bool isfinite() const
  {
    bool finite = true;
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < M; ++j)
        finite &= std::isfinite(eigen_()(i, j));
    return finite;
  }

  // Add another matrix.
  Matrix<N, M, T>& operator+=(Matrix<N, M, T> const& m2)
  {
    eigen_() += m2.eigen();
    return static_cast<Matrix<N, M, T>&>(*this);
  }

  // Subtract another matrix.
  Matrix<N, M, T>& operator-=(Matrix<N, M, T> const& m2)
  {
    eigen_() -= m2.eigen();
    return static_cast<Matrix<N, M, T>&>(*this);
  }

  // Multiply the matrix with a scalar.
  Matrix<N, M, T>& operator*=(T scalar)
  {
    eigen_() *= scalar;
    return static_cast<Matrix<N, M, T>&>(*this);
  }

  // Divide the matrix by a scalar.
  Matrix<N, M, T>& operator/=(T scalar)
  {
    eigen_() /= scalar;
    return static_cast<Matrix<N, M, T>&>(*this);
  }

  // Multiply with another matrix form the right.
  template<int P>
  Matrix<N, P, T> operator*(Matrix<M, P, T> const& m2) const
  {
    return {eigen_() * m2.eigen()};
  }

  // Multiply with a vector from the right.
  Vector<N, T> operator*(Vector<M, T> const& v) const
  {
    return {eigen_() * v.eigen()};
  }

  // Return the inverse of the matrix.
  Matrix<N, M, T> inverse() const requires (N == M) { return {eigen_().inverse()}; }

  // Return the transpose of the matrix.
  Matrix<M, N, T> transpose() const { return {eigen_().transpose()}; }

  // Negate this matrix.
  void negate() { eigen_() = -eigen_(); }

  // Return the negated matrix.
  Matrix<N, M, T> operator-() const { return {-eigen_()}; }

  // Add another matrix.
  Matrix<N, M, T> operator+(Matrix<N, M, T> const& m2) const { return {eigen_() + m2.eigen()}; }

  // Subtract another matrix.
  Matrix<N, M, T> operator-(Matrix<N, M, T> const& m2) const { return {eigen_() - m2.eigen()}; }

  // Multiply by a scalar.
  Matrix<N, M, T> operator*(T scalar) const { return {eigen_() * scalar}; }

  // Divide by a scalar.
  Matrix<N, M, T> operator/(T scalar) const { return {eigen_() / scalar}; }
};

// Matrix
//
// an N⨉M matrix
//
template<int N, int M, typename T = double>
class Matrix :
  public MatrixData<N, M, T>,               // Wraps underlying Eigen type and defines constructors.
  public MatrixOps<MatrixTypes<N, M, T>>    // Defines possible operations on the Matrix (uses the above specialization).
{
 public:
  using MatrixData<N, M, T>::MatrixData;
  using MatrixData<N, M, T>::operator=;
};

//-----------------------------------------------------------------------------
// Free function binary operators.
//

template<typename DerivedTypes>
DerivedTypes::template vector_type<DerivedTypes::cols>
operator*(typename DerivedTypes::template vector_type<DerivedTypes::rows> const& v, MatrixOps<DerivedTypes> const& m)
{
  return static_cast<DerivedTypes::template vector_type<DerivedTypes::cols>>(v.raw() * static_cast<DerivedTypes::derived_type const&>(m).raw());
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type operator*(typename DerivedTypes::scalar_type scalar, MatrixOps<DerivedTypes> const& m)
{
  return static_cast<DerivedTypes::derived_type>(scalar * static_cast<DerivedTypes::derived_type const&>(m).raw());
}

// Specialization of the free function binary operators.

template<int N, int M, typename T>
Vector<M, T> operator*(Vector<N, T> const& v, Matrix<N, M, T> const& m)
{
  return {v.eigen().transpose() * m.eigen()};
}

template<int N, int M, typename T>
Matrix<N, M, T> operator*(T scalar, MatrixOps<MatrixTypes<N, M, T>> const& m)
{
  return {scalar * static_cast<Matrix<N, M, T> const&>(m).eigen()};
}

#ifdef CWDEBUG
template<int N, int M, typename T>
void MatrixData<N, M, T>::print_on(std::ostream& os) const
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

} // namespace math

#endif // MATH_MATRIX_H
