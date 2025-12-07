#ifndef MATH_VECTOR_H
#define MATH_VECTOR_H

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <type_traits>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif
#if CW_DEBUG
#include "utils/almost_equal.h"
#include <limits>
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<int N, typename T>
class Point;

template<int N, typename T>
class Direction;

template<int N, typename T>
class LinePiece;

template<int N, typename T = double>
class Vector
{
 public:
  using scalar_type = T;
  using eigen_type = Eigen::Matrix<T, N, 1>;
  using point_type = Point<N, T>;
  using direction_type = Direction<N, T>;
  using line_piece_type = LinePiece<N, T>;

 protected:
  eigen_type v_;

 public:
  // Construct an uninitialized Vector.
  Vector() = default;

  // Construct a standard basis vector e_i.
  Vector(int i)
  {
    ASSERT(0 <= i && i < N);
    v_.setZero();
    v_[i] = T{1};
  }

  // Construct a vector from its x,y,... coordinates.
  template<typename... U>
    requires (sizeof...(U) == N &&                                              // Exactly N coefficients.
      (std::convertible_to<U, T> && ...) &&                                     // All convertible to T.
      !(N == 1 && ((std::same_as<std::remove_cvref_t<U>, Vector> || ...))))     // Do not replace copy/move constructor.
  constexpr Vector(U&&... xs) : v_(static_cast<T>(std::forward<U>(xs))...) { }

  // Construct a Vector that points in direction and has length.
  // Also used for automatic conversion from a Direction to a Vector.
  Vector(direction_type direction, T length = 1.0) : v_(direction.eigen() * length) { }

  // Construct a Vector from two points. If the second point is not given it defaults to the origin.
  // The direction is from the second argument (or origin) to the first argument.
  Vector(point_type const& from, point_type const& to) : v_(to.eigen() - from.eigen()) { }
  explicit Vector(point_type const& to) : v_(to.eigen()) { }

  // Construct a Vector from a Line, pointing from the first point to the second point.
  explicit Vector(line_piece_type const& line_piece) : Vector(line_piece.from(), line_piece.to()) { }

  // Construct a Vector from it eigen_type.
  Vector(eigen_type const& v) : v_(v) { }

  eigen_type& eigen() { return v_; }
  eigen_type const& eigen() const { return v_; }

  T& operator[](int i) { return v_(i); }
  T const& operator[](int i) const { return v_(i); }
  T& x() { static_assert(N >= 1); return v_[0]; }
  T& y() { static_assert(N >= 2); return v_[1]; }
  T& z() { static_assert(N >= 3); return v_[2]; }
  T x() const { static_assert(N >= 1); return v_[0]; }
  T y() const { static_assert(N >= 2); return v_[1]; }
  T z() const { static_assert(N >= 3); return v_[2]; }

  // Return dot product with v2.
  T dot(Vector const& v2) const { return v_.dot(v2.v_); }

  // Return the cross product with v2.
  auto cross(Vector const& v2) const requires (N == 2 || N == 3)
  {
    if constexpr (N == 2)
      return x() * v2.y() - y() * v2.x();
    else // N == 3
      return Vector{v_.cross(v2.v_)};
  }

  // Construct a Direction from this vector.
  direction_type direction() const { return direction_type{point_type{v_}}; }

  // Return the length of the vector.
  T norm() const { return v_.norm(); }

  // Return the square of the length of the vector.
  T norm_squared() const { return v_.squaredNorm(); }

  // Convert the vector to a point.
  point_type as_point() const { return {v_}; }

  bool isnan() const
  {
    bool has_nan = false;
    for (int i = 0; i < N; ++i)
      has_nan |= std::isnan(v_(i));
    return has_nan;
  }

  // Check if all coordinates are finite (not infinite or NaN).
  bool isfinite() const
  {
    bool finite = true;
    for (int i = 0; i < N; ++i)
      finite &= std::isfinite(v_(i));
    return finite;
  }

  // Returns true upon success.
  bool normalize()
  {
    // Turn the vector into a unit vector.
    v_.normalize();
    if (!isfinite())
      return false;
    ASSERT(utils::almost_equal(v_.squaredNorm(), 1.0, 1024.0 * std::numeric_limits<T>::epsilon()));
    return true;
  }

 public:
  // Return the vector rotated 90 degrees counter-clockwise.
  Vector rotate_90_degrees() const requires (N == 2) { return { -y(), x() }; }

  // Return the vector rotated 180 degrees.
  Vector rotate_180_degrees() const requires (N == 2)  { return { -x(), -y() }; }

  // Return the vector rotated 270 degrees.
  Vector rotate_270_degrees() const requires (N == 2) { return { y(), -x() }; }

  // Add another vector.
  Vector& operator+=(Vector const& v2)
  {
    v_ += v2.v_;
    return *this;
  }

  // Subtract another vector.
  Vector& operator-=(Vector const& v2)
  {
    v_ -= v2.v_;
    return *this;
  }

  // Multiply the vector with a scalar.
  Vector& operator*=(T scalar)
  {
    v_ *= scalar;
    return *this;
  }

  // Divide the vector by a scalar.
  Vector& operator/=(T scalar)
  {
    v_ /= scalar;
    return *this;
  }

  // Negate this vector.
  void negate()
  {
    v_ = -v_;
  }

  // Return the negated vector.
  Vector operator-() const
  {
    return {-v_};
  }

  // Divide by a scalar.
  Vector operator/(double scalar) const
  {
    return {v_ / scalar};
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << '[';
    char const* separator = "";
    for (int i = 0; i < N; ++i)
    {
      os << separator << v_[i];
      separator = ", ";
    }
    os << ']';
  }
#endif
};

template<int N, typename T>
Vector<N, T> operator*(double length, Vector<N, T> const& v2)
{
  return {length * v2.eigen()};
}

template<int N, typename T>
Point<N, T> operator+(Point<N, T> const& point, Vector<N, T> const& v2)
{
  return {point.eigen() + v2.eigen()};
}

template<int N, typename T>
Point<N, T>operator-(Point<N, T> const& point, Vector<N, T> const& v2)
{
  return {point.eigen() - v2.eigen()};
}

template<int N, typename T>
Vector<N, T> operator+(Vector<N, T> const& v1, Vector<N, T> const& v2)
{
  return {v1.eigen() + v2.eigen()};
}

template<int N, typename T>
Vector<N, T> operator-(Vector<N, T> const& v1, Vector<N, T> const& v2)
{
  return {v1.eigen() - v2.eigen()};
}

} // namespace math

#endif // MATH_VECTOR_H
