#ifndef MATH_DIRECTION_H
#define MATH_DIRECTION_H

#include <cmath>
#include <limits>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif
#include "debug.h"

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<int N, typename T>
class LinePiece;

template<int N, typename T>
class Line;

template<int N, typename T>
class Direction;

template<typename DerivedTypes>
struct DirectionOps;

template<int N, typename T>
class DirectionData
{
 public:
  using eigen_type = Eigen::Matrix<T, N, 1>;

 protected:
  eigen_type d_;                // d_ is a unit vector pointing in the intended direction.

 public:
  // Construct an uninitialized Direction.
  DirectionData() = default;

  // Construct a Direction that points in the direction theta (in radians): the counter-clockwise angle with the positive x-axis.
  DirectionData(T theta) requires (N == 2) : d_(std::cos(theta), std::sin(theta)) { }

  // Construct a Direction from two points.
  explicit DirectionData(Point<N, T> const& from, Point<N, T> const& to) : d_(to.eigen() - from.eigen())
  {
    // Normalize the unit vector.
    d_.normalize();
#if CW_DEBUG
    bool has_nan = false;
    for (int i = 0; i < N; ++i)
      has_nan |= std::isnan(d_(i));
    ASSERT(!has_nan);
#endif
  }

  // If only one point is give, the direction is from the origin to that point.
  explicit DirectionData(Point<N, T> const& to) : DirectionData(Point<N, T>(Point<N, T>::eigen_type::Zero()), to) { }

  // Construct a Direction from a LinePiece.
  DirectionData(LinePiece<N, T> const& line_piece);

  // Construct a Direction from a Line.
  DirectionData(Line<N, T> const& line);

  // Accessors.
  eigen_type& eigen() { return d_; }
  eigen_type const& eigen() const { return d_; }

 protected:
  // For normal() and inverse().
  template<typename DerivedTypes>
  friend struct DirectionOps;
  template<typename... U>
    requires (sizeof...(U) == N &&                                                      // Exactly N coefficients.
      (std::convertible_to<U, T> && ...) &&                                             // All convertible to T.
      !(N == 1 && ((std::same_as<std::remove_cvref_t<U>, DirectionData> || ...))))      // Do not replace copy/move constructor.
  constexpr DirectionData(U&&... xs) : d_(static_cast<T>(std::forward<U>(xs))...) { }

  DirectionData(eigen_type d) : d_(d) { }

#ifdef CWDEBUG
 public:
  void print_on(std::ostream& os) const;
#endif
};

template<typename DerivedTypes>
struct DirectionOps
{
 private:
  static constexpr int derived_n = DerivedTypes::n;
  using derived_scalar_type    = typename DerivedTypes::scalar_type;
  using derived_type           = typename DerivedTypes::derived_type;

  auto& raw_() { return static_cast<derived_type*>(this)->raw(); }
  auto const& raw_() const { return static_cast<derived_type const*>(this)->raw(); }

 public:
  derived_scalar_type const& operator[](int i) const { return raw_().operator[](i); }
  derived_scalar_type x() const { return raw_().x(); }
  derived_scalar_type y() const { return raw_().y(); }
  derived_scalar_type z() const { return raw_().z(); }

  inline derived_scalar_type dot(derived_type const& d2) const;
  inline derived_scalar_type as_angle(int x = 0, int y = 1) const;

  inline derived_type normal() const;
  inline derived_type inverse() const;
  inline derived_type normal_inverse() const;

  static derived_type const& up;
  static derived_type const& down;
  static derived_type const& left;
  static derived_type const& right;
};

//static
template<typename DerivedTypes>
typename DerivedTypes::derived_type const& DirectionOps<DerivedTypes>::up = reinterpret_cast<derived_type const&>(Direction<derived_n, derived_scalar_type>::up);

//static
template<typename DerivedTypes>
typename DerivedTypes::derived_type const& DirectionOps<DerivedTypes>::down = reinterpret_cast<derived_type const&>(Direction<derived_n, derived_scalar_type>::down);

//static
template<typename DerivedTypes>
typename DerivedTypes::derived_type const& DirectionOps<DerivedTypes>::left = reinterpret_cast<derived_type const&>(Direction<derived_n, derived_scalar_type>::left);

//static
template<typename DerivedTypes>
typename DerivedTypes::derived_type const& DirectionOps<DerivedTypes>::right = reinterpret_cast<derived_type const&>(Direction<derived_n, derived_scalar_type>::right);

// Forward declaration, required for DirectionTypes<N, T>::derived_type.
template<int N, typename T>
class Direction;

template<int N, typename T>
struct DirectionTypes
{
  static constexpr int n = N;
  using scalar_type    = T;
  using derived_type   = Direction<N, T>;
};

// Specialization of the Direction operators specifically for math::Direction itself.
template<int N, typename T>
struct DirectionOps<DirectionTypes<N, T>>
{
 private:
  // Get underlying Eigen type.
  auto& eigen_() { return static_cast<Direction<N, T>*>(this)->eigen(); }
  auto const& eigen_() const { return static_cast<Direction<N, T> const*>(this)->eigen(); }

 public:
  T const& operator[](int i) const { return eigen_()(i); }
  T x() const { static_assert(N >= 1); return eigen_()[0]; }
  T y() const { static_assert(N >= 2); return eigen_()[1]; }
  T z() const { static_assert(N >= 3); return eigen_()[2]; }

  // Return dot product with d2.
  inline T dot(Direction<N, T> const& d2) const;

  // Returns the polar angle of the projection onto the (x,y) plane, in (-π, π] radians.
  inline T as_angle(int x = 0, int y = 1) const requires (N >= 2);

  // Return the direction rotated 90 degrees counter-clockwise.
  inline Direction<N, T> normal() const requires (N == 2);

  // Return the inverse of the direction (aka rotated 180 degrees if N=2).
  inline Direction<N, T> inverse() const;

  // Return the direction rotated 270 degrees.
  inline Direction<N, T> normal_inverse() const requires (N == 2);

  // A few convenience directions.
  static Direction<N, T> const up;
  static Direction<N, T> const down;
  static Direction<N, T> const left;
  static Direction<N, T> const right;
};

template<int N, typename T = double>
class Direction :
  public DirectionData<N, T>,                   // Wraps underlying Eigen type and defines constructors.
  public DirectionOps<DirectionTypes<N, T>>     // Defines possible operations on the Direction (uses the above specialization).
{
 public:
  using DirectionData<N, T>::DirectionData;
};

} // namespace math

#endif // MATH_DIRECTION_H

#ifndef MATH_LINE_PIECE_H
#include "LinePiece.h"
#endif // MATH_LINE_PIECE_H

#ifndef MATH_LINE_H
#include "Line.h"
#endif // MATH_LINE_H

#ifndef MATH_DIRECTION_H_definitions
#define MATH_DIRECTION_H_definitions

namespace math {

template<typename DerivedTypes>
typename DerivedTypes::scalar_type DirectionOps<DerivedTypes>::dot(derived_type const& d2) const
{
  return raw_().dot(d2.raw());
}

template<typename DerivedTypes>
typename DerivedTypes::scalar_type DirectionOps<DerivedTypes>::as_angle(int x, int y) const
{
  return raw_().as_angle(x, y);
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type DirectionOps<DerivedTypes>::normal() const
{
  return static_cast<derived_type>(raw_().normal());
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type DirectionOps<DerivedTypes>::inverse() const
{
  return static_cast<derived_type>(raw_().inverse());
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type DirectionOps<DerivedTypes>::normal_inverse() const
{
  return static_cast<derived_type>(raw_().normal_inverse());
}

// Specialization of the Direction operators specifically for math::Direction itself.

template<int N, typename T>
T DirectionOps<DirectionTypes<N, T>>::dot(Direction<N, T> const& d2) const
{
  return eigen_().dot(d2.eigen());
}

template<int N, typename T>
T DirectionOps<DirectionTypes<N, T>>::as_angle(int x, int y) const requires (N >= 2)
{
  ASSERT(0 <= x && x < N);
  ASSERT(0 <= y && y < N);
  ASSERT(x != y);

  T const px = eigen_()(x);
  T const py = eigen_()(y);

  // Degenerate if projection is zero-length.
  if (std::hypot(px, py) == T(0)) [[unlikely]]
    return std::numeric_limits<T>::quiet_NaN();

  return std::atan2(py, px);
}

template<int N, typename T>
Direction<N, T> DirectionOps<DirectionTypes<N, T>>::normal() const requires (N == 2)
{
  return { -y(), x() };
}

template<int N, typename T>
Direction<N, T> DirectionOps<DirectionTypes<N, T>>::inverse() const
{
  return {-eigen_()};
}

template<int N, typename T>
Direction<N, T> DirectionOps<DirectionTypes<N, T>>::normal_inverse() const requires (N == 2)
{
  return { y(), -x() };
}

//static
template<int N, typename T>
Direction<N, T> const DirectionOps<DirectionTypes<N, T>>::up{0, 1};
//static
template<int N, typename T>
Direction<N, T> const DirectionOps<DirectionTypes<N, T>>::down{0, -1};
//static
template<int N, typename T>
Direction<N, T> const DirectionOps<DirectionTypes<N, T>>::left{-1, 0};
//static
template<int N, typename T>
Direction<N, T> const DirectionOps<DirectionTypes<N, T>>::right{1, 0};

template<int N, typename T>
DirectionData<N, T>::DirectionData(LinePiece<N, T> const& line_piece) : DirectionData(line_piece.from(), line_piece.to())
{
}

template<int N, typename T>
DirectionData<N, T>::DirectionData(Line<N, T> const& line) : DirectionData(line.direction())
{
}

#ifdef CWDEBUG
template<int N, typename T>
void DirectionData<N, T>::print_on(std::ostream& os) const
{
  os << "{d_:" << d_ << '}';
}
#endif

} // namespace math

#endif // MATH_DIRECTION_H_definitions
