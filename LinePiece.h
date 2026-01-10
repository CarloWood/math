#ifndef MATH_LINE_PIECE_H
#define MATH_LINE_PIECE_H

#include "Point.h"
#include <cmath>

namespace math {

template<int N, typename T>
class Direction;

// Friend of LinePieceData.
template<typename DerivedTypes>
struct LinePieceOps;

template<int N, typename T>
class LinePieceData
{
 public:
  using point_type = Point<N, T>;

 protected:
  point_type from_;
  point_type to_;

 protected:
  template<typename DerivedTypes>
  friend struct LinePieceOps;

 public:
  // Construct an uninitialized line.
  LinePieceData() = default;

  // Construct a LinePiece from two points.
  LinePieceData(point_type const& from, point_type const& to) : from_(from), to_(to) { }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

template<typename DerivedTypes>
struct LinePieceOps
{
 private:
  using derived_scalar_type    = typename DerivedTypes::scalar_type;
  using derived_type           = typename DerivedTypes::derived_type;
  using derived_point_type     = typename DerivedTypes::point_type;
  using derived_direction_type = typename DerivedTypes::direction_type;

  auto& raw_() { return static_cast<derived_type*>(this)->raw(); }
  auto const& raw_() const { return static_cast<derived_type const*>(this)->raw(); }

 public:
  derived_point_type from() const { return static_cast<derived_point_type>(raw_().from()); }
  derived_point_type to() const { return static_cast<derived_point_type>(raw_().to()); }
  derived_scalar_type norm() const { return raw_().norm(); }
  derived_scalar_type norm_squared() const { return raw_().norm_squared(); }
  derived_direction_type direction() const { return static_cast<derived_direction_type>(raw_().direction()); }
};

template<int N, typename T>
struct LinePieceTypes
{
  using scalar_type    = T;
  using derived_type   = LinePiece<N, T>;
  using point_type     = Point<N, T>;
  using direction_type = Direction<N, T>;
};

// Specialization of LinePiece operations specifically for math::LinePiece itself.
template<int N, typename T>
struct LinePieceOps<LinePieceTypes<N, T>>
{
 public:
  // Accessors.
  Point<N, T> const& from() const { return static_cast<LinePiece<N, T> const*>(this)->from_; }
  Point<N, T> const& to() const { return static_cast<LinePiece<N, T> const*>(this)->to_; }

  T norm() const { return (to().eigen() - from().eigen()).norm(); }
  T norm_squared() const { return (to().eigen() - from().eigen()).squaredNorm(); }

  Direction<N, T> direction() const;
};

template<int N, typename T = double>
class LinePiece :
  public LinePieceData<N, T>,                    // Wraps stored endpoints and defines constructors.
  public LinePieceOps<LinePieceTypes<N, T>>      // Defines possible operations on the LinePiece (uses the above specialization).
{
 public:
  using LinePieceData<N, T>::LinePieceData;
};

// Implicit class template argument deduction guides are not generated from inherited constructors (LinePieceData),
// therefore we have to add explicit deduction guides back for LinePiece.
template<int N, class T>
LinePiece(Point<N, T> const&, Point<N, T> const&) -> LinePiece<N, T>;

} // namespace math

#endif // MATH_LINE_PIECE_H

#ifndef MATH_DIRECTION_H
#include "Direction.h"
#endif // MATH_DIRECTION_H

#ifndef MATH_LINE_PIECE_H_definitions
#define MATH_LINE_PIECE_H_definitions

namespace math {

template<int N, typename T>
Direction<N, T> LinePieceOps<LinePieceTypes<N, T>>::direction() const
{
  return Direction<N, T>{from(), to()};
}

#ifdef CWDEBUG
template<int N, typename T>
void LinePieceData<N, T>::print_on(std::ostream& os) const
{
  os << "{from:" << from_ << ", to:" << to_ << "}";
}
#endif

} // namespace math

#endif // MATH_LINE_PIECE_H_definitions
