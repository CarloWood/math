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

template<int N, typename T>
class Vector;

template<typename DerivedTypes>
struct VectorOps;

template<int N, typename T>
class VectorData
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
  VectorData() = default;

  // Construct a standard basis vector e_i, or zero if i == -1.
  VectorData(int i)
  {
    ASSERT(-1 <= i && i < N);
    v_.setZero();
    if (i != -1)
      v_[i] = T{1};
  }

  // Construct a vector from its x,y,... coordinates.
  template<typename... U>
    requires (sizeof...(U) == N &&                                              // Exactly N coefficients.
      (std::convertible_to<U, T> && ...) &&                                     // All convertible to T.
      !(N == 1 && ((std::same_as<std::remove_cvref_t<U>, VectorData> || ...)))) // Do not replace copy/move constructor.
  constexpr VectorData(U&&... xs) : v_(static_cast<T>(std::forward<U>(xs))...) { }

  // Construct a Vector that points in direction and has length.
  // Also used for automatic conversion from a Direction to a Vector.
  VectorData(direction_type direction, T length = 1.0) : v_(direction.eigen() * length) { }

  // Construct a Vector from two points. If the second point is not given it defaults to the origin.
  // The direction is from the second argument (or origin) to the first argument.
  VectorData(point_type const& from, point_type const& to) : v_(to.eigen() - from.eigen()) { }
  explicit VectorData(point_type const& to) : v_(to.eigen()) { }

  // Construct a Vector from a Line, pointing from the first point to the second point.
  explicit VectorData(line_piece_type const& line_piece) : VectorData(line_piece.from(), line_piece.to()) { }

  // Construct a Vector from it eigen_type.
  VectorData(eigen_type const& v) : v_(v) { }

  eigen_type& eigen() { return v_; }
  eigen_type const& eigen() const { return v_; }

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

template<typename DerivedTypes>
struct VectorOps
{
 private:
  static constexpr int derived_n = DerivedTypes::n;
  using derived_scalar_type    = typename DerivedTypes::scalar_type;
  using derived_type           = typename DerivedTypes::derived_type;
  using derived_point_type     = typename DerivedTypes::point_type;
  using derived_direction_type = typename DerivedTypes::direction_type;

  auto& raw_() { return static_cast<derived_type*>(this)->raw(); }
  auto const& raw_() const { return static_cast<derived_type const*>(this)->raw(); }

 public:
  derived_scalar_type&       operator[](int i)       { return raw_().operator[](i); }
  derived_scalar_type const& operator[](int i) const { return raw_().operator[](i); }
  derived_scalar_type& x()      { return raw_().x(); }
  derived_scalar_type& y()      { return raw_().y(); }
  derived_scalar_type& z()      { return raw_().z(); }
  derived_scalar_type x() const { return raw_().x(); }
  derived_scalar_type y() const { return raw_().y(); }
  derived_scalar_type z() const { return raw_().z(); }

  derived_scalar_type dot(derived_type const& v2) const { return raw_().dot(v2.raw()); }
  auto cross(derived_type const& v2) const requires (derived_n == 2 || derived_n == 3)
  {
    if constexpr (derived_n == 2)
      return raw_().cross(v2.raw());
    else
      return static_cast<derived_type>(raw_().cross(v2.raw()));
  }

  derived_direction_type direction() const { return static_cast<derived_direction_type>(raw_().direction()); }
  derived_scalar_type norm() const { return raw_().norm(); }
  derived_scalar_type norm_squared() const { return raw_().norm_squared(); }
  derived_point_type as_point() const { return static_cast<derived_point_type>(raw_().as_point()); }
  bool isnan() const { return raw_().isnan(); }
  bool isfinite() const { return raw_().isfinite(); }
  bool normalize() { return raw_().normalize(); }

  derived_type rotate_90_degrees() const requires (derived_n == 2) { return static_cast<derived_type>(raw_().rotate_90_degrees()); }
  derived_type rotate_180_degrees() const requires (derived_n == 2) { return static_cast<derived_type>(raw_().rotate_180_degrees()); }
  derived_type rotate_270_degrees() const requires (derived_n == 2) { return static_cast<derived_type>(raw_().rotate_270_degrees()); }

  derived_type& operator+=(derived_type const& v2) { raw_().operator+=(v2.raw()); return static_cast<derived_type&>(*this); }
  derived_type& operator-=(derived_type const& v2) { raw_().operator-=(v2.raw()); return static_cast<derived_type&>(*this); }
  derived_type& operator*=(derived_scalar_type scalar) { raw_().operator*=(scalar); return static_cast<derived_type&>(*this); }
  derived_type& operator/=(derived_scalar_type scalar) { raw_().operator/=(scalar); return static_cast<derived_type&>(*this); }

  void negate() { raw_().negate(); }
  derived_type operator-() const { return static_cast<derived_type>(raw_().operator-()); }
  derived_type operator/(double scalar) const { return static_cast<derived_type>(raw_().operator/(scalar)); }

  friend derived_type operator*(double length, VectorOps const& v2) { return derived_type{length * v2.raw_()}; }
  friend derived_type operator+(VectorOps const& v1, VectorOps const& v2) { return derived_type{v1.raw_() + v2.raw_()}; }
  friend derived_type operator-(VectorOps const& v1, VectorOps const& v2) { return derived_type{v1.raw_() - v2.raw_()}; }
};

template<int N, typename T>
struct VectorTypes
{
  static constexpr int n = N;
  using scalar_type    = T;
  using point_type     = Point<N, T>;
  using direction_type = Direction<N, T>;
  using derived_type   = Vector<N, T>;
};

// Specialization of the Vector operators specifically for math::Vector itself.
template<int N, typename T>
struct VectorOps<VectorTypes<N, T>>
{
 private:
  // Get underlying Eigen type.
  auto& eigen_() { return static_cast<Vector<N, T>*>(this)->eigen(); }
  auto const& eigen_() const { return static_cast<Vector<N, T> const*>(this)->eigen(); }

 public:
  T&       operator[](int i)       { return eigen_()(i); }
  T const& operator[](int i) const { return eigen_()(i); }
  T& x() { static_assert(N >= 1); return eigen_()[0]; }
  T& y() { static_assert(N >= 2); return eigen_()[1]; }
  T& z() { static_assert(N >= 3); return eigen_()[2]; }
  T x() const { static_assert(N >= 1); return eigen_()[0]; }
  T y() const { static_assert(N >= 2); return eigen_()[1]; }
  T z() const { static_assert(N >= 3); return eigen_()[2]; }

  // Return dot product with v2.
  T dot(Vector<N, T> const& v2) const { return eigen_().dot(v2.eigen()); }

  // Return the cross product with v2.
  auto cross(Vector<N, T> const& v2) const requires (N == 2 || N == 3)
  {
    if constexpr (N == 2)
      return x() * v2.y() - y() * v2.x();
    else // N == 3
      return Vector<N, T>{eigen_().cross(v2.eigen())};
  }

  // Construct a Direction from this vector.
  Direction<N, T> direction() const { return Direction<N, T>{Point<N, T>{eigen_()}}; }

  // Return the length of the vector.
  T norm() const { return eigen_().norm(); }

  // Return the square of the length of the vector.
  T norm_squared() const { return eigen_().squaredNorm(); }

  // Convert the vector to a point.
  Point<N, T> as_point() const { return {eigen_()}; }

  bool isnan() const
  {
    bool has_nan = false;
    for (int i = 0; i < N; ++i)
      has_nan |= std::isnan(eigen_()(i));
    return has_nan;
  }

  // Check if all coordinates are finite (not infinite or NaN).
  bool isfinite() const
  {
    bool finite = true;
    for (int i = 0; i < N; ++i)
      finite &= std::isfinite(eigen_()(i));
    return finite;
  }

  // Returns true upon success.
  bool normalize()
  {
    // Turn the vector into a unit vector.
    eigen_().normalize();
    if (!isfinite())
      return false;
    ASSERT(utils::almost_equal(eigen_().squaredNorm(), 1.0, 1024.0 * std::numeric_limits<T>::epsilon()));
    return true;
  }

  // Return the vector rotated 90 degrees counter-clockwise.
  Vector<N, T> rotate_90_degrees() const requires (N == 2) { return { -y(), x() }; }

  // Return the vector rotated 180 degrees.
  Vector<N, T> rotate_180_degrees() const requires (N == 2) { return { -x(), -y() }; }

  // Return the vector rotated 270 degrees.
  Vector<N, T> rotate_270_degrees() const requires (N == 2) { return { y(), -x() }; }

  // Add another vector.
  Vector<N, T>& operator+=(Vector<N, T> const& v2)
  {
    eigen_() += v2.eigen();
    return static_cast<Vector<N, T>&>(*this);
  }

  // Subtract another vector.
  Vector<N, T>& operator-=(Vector<N, T> const& v2)
  {
    eigen_() -= v2.eigen();
    return static_cast<Vector<N, T>&>(*this);
  }

  // Multiply the vector with a scalar.
  Vector<N, T>& operator*=(T scalar)
  {
    eigen_() *= scalar;
    return static_cast<Vector<N, T>&>(*this);
  }

  // Divide the vector by a scalar.
  Vector<N, T>& operator/=(T scalar)
  {
    eigen_() /= scalar;
    return static_cast<Vector<N, T>&>(*this);
  }

  // Negate this vector.
  void negate() { eigen_() = -eigen_(); }

  // Return the negated vector.
  Vector<N, T> operator-() const { return {-eigen_()}; }

  // Divide by a scalar.
  Vector<N, T> operator/(double scalar) const { return {eigen_() / scalar}; }

  // Binary operators.
  friend Vector<N, T> operator*(double length, VectorOps const& v2) { return {length * v2.eigen_()}; }
  friend Vector<N, T> operator+(VectorOps const& v1, VectorOps const& v2) { return {v1.eigen_() + v2.eigen_()}; }
  friend Vector<N, T> operator-(VectorOps const& v1, VectorOps const& v2) { return {v1.eigen_() - v2.eigen_()}; }
};

template<int N, typename T = double>
class Vector :
  public VectorData<N, T>,               // Wraps underlying Eigen type and defines constructors.
  public VectorOps<VectorTypes<N, T>>    // Defines possible operations on the Vector (uses the above specialization).
{
 public:
  using VectorData<N, T>::VectorData;
};

// Implicit class template argument deduction guides are not generated from inherited constructors (VectorData),
// therefore we have to add explicit deduction guides back for Vector where possible.
template<int N, class T>
Vector(Direction<N, T>, T) -> Vector<N, T>;

template<int N, class T>
Vector(Direction<N, T>) -> Vector<N, T>;

template<int N, class T>
Vector(Point<N, T> const&, Point<N, T> const&) -> Vector<N, T>;

template<int N, class T>
Vector(Point<N, T> const&) -> Vector<N, T>;

template<int N, class T>
Vector(LinePiece<N, T> const&) -> Vector<N, T>;

} // namespace math

#endif // MATH_VECTOR_H
