#ifndef MATH_LINE_PIECE_H
#define MATH_LINE_PIECE_H

#include "Point.h"
#include <cmath>

namespace math {

template<int N, typename T>
class Direction;

template<int N, typename T = double>
class LinePiece
{
 public:
  using scalar_type = T;
  using point_type = Point<N, T>;
  using direction_type = Direction<N, T>;

 protected:
  point_type from_;
  point_type to_;

 public:
  LinePiece() = default;
  LinePiece(point_type const& from, point_type const& to) : from_(from), to_(to) { }

  point_type const& from() const { return from_; }
  point_type const& to() const { return to_; }
  T norm() const { return (to_.eigen() - from_.eigen()).norm(); }

  direction_type direction() const;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

} // namespace math

#endif // MATH_LINE_PIECE_H

#ifndef MATH_DIRECTION_H
#include "Direction.h"
#endif // MATH_DIRECTION_H

#ifndef MATH_LINE_PIECE_H_definitions
#define MATH_LINE_PIECE_H_definitions

namespace math {

template<int N, typename T>
Direction<N, T> LinePiece<N, T>::direction() const
{
  return Direction{from_, to_};
}

#ifdef CWDEBUG
template<int N, typename T>
void LinePiece<N, T>::print_on(std::ostream& os) const
{
  os << "{from:" << from_ << ", to:" << to_ << "}";
}
#endif

} // namespace math

#endif // MATH_LINE_PIECE_H_definitions
