#ifndef MATH_LINE_H
#define MATH_LINE_H

#include "Direction.h"
#include "Point.h"

namespace math {

template<int N, typename T>
class LineData
{
 public:
  using point_type = Point<N, T>;
  using direction_type = Direction<N, T>;

 protected:
  point_type point_;
  direction_type direction_;

 public:
  // Construct an uninitialized line.
  LineData() = default;

  // Construct a line through point with direction.
  LineData(point_type const& point, direction_type const& direction) : point_(point), direction_(direction) { }

  point_type const& point() const { return point_; }
  direction_type const& direction() const { return direction_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

template<typename DerivedTypes>
struct LineOps
{
 private:
  static constexpr int derived_n = DerivedTypes::n;
  using derived_type           = typename DerivedTypes::derived_type;
  using derived_point_type     = typename DerivedTypes::point_type;
  using derived_direction_type = typename DerivedTypes::direction_type;

  auto& raw_() { return static_cast<derived_type*>(this)->raw(); }
  auto const& raw_() const { return static_cast<derived_type const*>(this)->raw(); }

 public:
  derived_point_type const& point() const { return raw_().point(); }
  derived_direction_type const& direction() const { return raw_().direction(); }

  operator derived_direction_type const&() const { return direction(); }

  derived_point_type intersection_with(derived_type const& line2) const requires (derived_n == 2)
  {
    return static_cast<derived_point_type>(raw_().intersection_with(line2.raw()));
  }
};

// Forward declaration, required for LineTypes<N, T>::derived_type.
template<int N, typename T>
class Line;

template<int N, typename T>
struct LineTypes
{
  static constexpr int n = N;
  using scalar_type    = T;
  using derived_type   = Line<N, T>;
  using point_type     = Point<N, T>;
  using direction_type = Direction<N, T>;
};

// Specialization of Line operations specifically for math::Line itself.
template<int N, typename T>
struct LineOps<LineTypes<N, T>>
{
 private:
  // Get underlying point and direction types.
  auto const& point_() const { return static_cast<Line<N, T> const*>(this)->point(); }
  auto const& direction_() const { return static_cast<Line<N, T> const*>(this)->direction(); }

 public:
  // It is allowed to pass a Line to a function that takes a Direction.
  operator Direction<N, T> const&() const { return direction_(); }

  // Return the intersection with another Line.
  Point<N, T> intersection_with(Line<N, T> const& line2) const requires (N == 2);
};

template<int N, typename T = double>
class Line :
  public LineData<N, T>,               // Wraps stored point and direction and defines constructors.
  public LineOps<LineTypes<N, T>>      // Defines possible operations on the Line (uses the above specialization).
{
 public:
  using LineData<N, T>::LineData;
};

// Implicit class template argument deduction guides are not generated from inherited constructors (LineData),
// therefore we have to add explicit deduction guides back for Line.
template<int N, class T>
Line(Point<N, T> const&, Direction<N, T> const&) -> Line<N, T>;

template<int N, typename T>
Point<N, T> LineOps<LineTypes<N, T>>::intersection_with(Line<N, T> const& L1) const requires (N == 2)
{
  // Line0: P0 + λ D0
  // Line1: P1 + ξ D1

  Point<N, T> const& P0 = point_();
  Direction<N, T> const& D0 = direction_();
  Point<N, T> const& P1 = L1.point();
  Direction<N, T> const& D1 = L1.direction();

  // Let N1 be D1 rotated counter-clockwise by PI/2 (this is floating-point round off error free).
  Direction<N, T> N1 = D1.normal();

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
  T P1P0_dot_N1 = (L1.point().x() - point_().x()) * N1.x() + (L1.point().y() - point_().y()) * N1.y();

  // Calculate lambda.
  T lambda = P1P0_dot_N1 / D0_dot_N1;

  // Return intersection point.
  return {point_().x() + lambda * direction_().x(), point_().y() + lambda * direction_().y()};
}

#ifdef CWDEBUG
template<int N, typename T>
void LineData<N, T>::print_on(std::ostream& os) const
{
  os << "{point:" << point_ << ", direction:" << direction_ << "}";
}
#endif

} // namespace math

#endif // MATH_LINE_H
