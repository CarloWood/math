#ifndef MATH_LINE_H
#define MATH_LINE_H

#include "Direction.h"
#include "Point.h"

namespace math {

template<int N, typename T = double>
class Line
{
 public:
  using point_type = Point<N, T>;
  using direction_type = Direction<N, T>;

 protected:
  point_type point_;
  direction_type direction_;

 public:
  // Construct an undefined line.
  Line() = default;
  // Construct a line through point with direction.
  Line(point_type const& point, direction_type const& direction) : point_(point), direction_(direction) { }

  point_type const& point() const { return point_; }
  direction_type const& direction() const { return direction_; }

  operator direction_type const&() const { return direction_; }

  point_type intersection_with(Line const& line2) const requires (N == 2);

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

template<int N, typename T>
Point<N, T> Line<N, T>::intersection_with(Line const& L1) const requires (N == 2)
{
  // Line0: P0 + λ D0
  // Line1: P1 + ξ D1

  point_type const& P0 = point_;
  direction_type const& D0 = direction_;
  point_type const& P1 = L1.point();
  direction_type const& D1 = L1.direction();

  // Let N1 be D1 rotated counter-clockwise by PI/2 (this is floating-point round off error free).
  direction_type N1 = D1.normal();

  // Take dot product of D0 with N1:
  T D0_dot_N1 = D0.dot(N1);

  //            intersection
  //                 \ /         P1
  //  --------+-------+-<--------+
  //       ^  |      /    D1  1  |
  //       |  |     /λ           |N1
  //     a |  |_  _/            1|
  //       | ^|   /|             |
  //       |b|| 1/               v
  //       | || /D0
  //       v v|/
  //          +P0
  //         /
  // intersection: P0 + λ D0 = P1 + ξ D1 -->
  //
  //   P1 - P0 = λ D0 - ξ D1
  //
  // Take dot product with N1 (D1 rotated PI/2 radians):
  //
  // λ = (P1 - P0)·N1 / (D0·N1)
  //
  // Note, in the above picture: a = -(P1 - P0)·N1, and b = -D0·N1

  // Take dot product of P1-P0 with N1:
  T P1P0_dot_N1 = (L1.point().x() - point_.x()) * N1.x() + (L1.point().y() - point_.y()) * N1.y();

  // Calculate lambda.
  T lambda = P1P0_dot_N1 / D0_dot_N1;

  // Return intersection point.
  return {point_.x() + lambda * direction_.x(), point_.y() + lambda * direction_.y()};
}

#ifdef CWDEBUG
template<int N, typename T>
void Line<N, T>::print_on(std::ostream& os) const
{
  os << "{point:" << point_ << ", direction:" << direction_ << "}";
}
#endif

} // namespace math

#endif // MATH_LINE_H
