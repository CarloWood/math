#ifndef MATH_DIRECTION_H
#define MATH_DIRECTION_H

#include "LinePiece.h"
#include <cmath>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif
#include "debug.h"

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<int N, typename T>
class Line;

template<int N, typename T = double>
class Direction
{
 public:
  using scalar_type = T;
  using eigen_type = Eigen::Matrix<T, N, 1>;
  using point_type = Point<N, T>;
  using line_piece_type = LinePiece<N, T>;
  using line_type = Line<N, T>;

 protected:
  eigen_type d_;                // d_ is a unit vector pointing in the intended direction.

 public:
  // Construct an undefined Direction.
  Direction() = default;

  // Construct a Direction that points in the direction theta (in radians): an angle with the positive x-axis.
  Direction(T theta) requires (N == 2) : d_(std::cos(theta), std::sin(theta)) { }

  // Construct a Direction from two points. If the second point is not given it defaults to the origin.
  explicit Direction(point_type const& from, point_type const& to) : d_(to.eigen() - from.eigen())
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
  explicit Direction(point_type const& to) : Direction(point_type(point_type::eigen_type::Zero()), to) { }

  // Construct a Direction from a LinePiece, pointing from the first point to the second point.
  Direction(line_piece_type const& line_piece) : Direction(line_piece.from(), line_piece.to()) { }

  // Construct a Direction from a Line.
  Direction(line_type const& line);

  eigen_type& eigen() { return d_; }
  eigen_type const& eigen() const { return d_; }

  T const& operator[](int i) const { return d_(i); }
  T x() const { static_assert(N >= 1); return d_[0]; }
  T y() const { static_assert(N >= 2); return d_[1]; }
  T z() const { static_assert(N >= 3); return d_[2]; }

  // Return dot product with d2.
  T dot(Direction const& d2) const { return d_.dot(d2.d_); }

  // Returns the polar angle of the projection onto the (x,y) plane, in (-π, π] radians.
  T as_angle(int x = 0, int y = 1) const requires (N >= 2)
  {
    ASSERT(0 <= x && x < N);
    ASSERT(0 <= y && y < N);
    ASSERT(x != y);

    T const px = d_(x);
    T const py = d_(y);

    // Degenerate if projection is zero-length.
    if (std::hypot(px, py) == T(0)) [[unlikely]]
      return std::numeric_limits<T>::quiet_NaN();

    return std::atan2(py, px);
  }

 protected:
  // For normal() and inverse().
  template<typename... U>
    requires (sizeof...(U) == N &&                                              // Exactly N coefficients.
      (std::convertible_to<U, T> && ...) &&                                     // All convertible to T.
      !(N == 1 && ((std::same_as<std::remove_cvref_t<U>, Direction> || ...))))  // Do not replace copy/move constructor.
  constexpr Direction(U&&... xs) : d_(static_cast<T>(std::forward<U>(xs))...) { }

 public:
  // Return the direction rotated 90 degrees counter-clockwise.
  Direction normal() const requires (N == 2) { return { -y(), x() }; }

  // Return the direction rotated 180 degrees.
  Direction inverse() const requires (N == 2) { return { -x(), -y() }; }

  // Return the direction rotated 270 degrees.
  Direction normal_inverse() const requires (N == 2) { return { y(), -x() }; }

  static Direction const up;
  static Direction const down;
  static Direction const left;
  static Direction const right;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

} // namespace math

#endif // MATH_DIRECTION_H

#ifndef MATH_LINE_H
#include "Line.h"
#endif // MATH_LINE_H

#ifndef MATH_DIRECTION_H_definitions
#define MATH_DIRECTION_H_definitions

namespace math {

template<int N, typename T>
Direction<N, T>::Direction(line_type const& line) : Direction(line.direction())
{
}

//static
template<int N, typename T>
Direction<N, T> const Direction<N, T>::up{0, 1};
//static
template<int N, typename T>
Direction<N, T> const Direction<N, T>::down{0, -1};
//static
template<int N, typename T>
Direction<N, T> const Direction<N, T>::left{-1, 0};
//static
template<int N, typename T>
Direction<N, T> const Direction<N, T>::right{1, 0};

#ifdef CWDEBUG
template<int N, typename T>
void Direction<N, T>::print_on(std::ostream& os) const
{
  os << "{d_:" << d_ << '}';
}
#endif

} // namespace math

#endif // MATH_DIRECTION_H_definitions
