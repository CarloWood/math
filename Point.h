#ifndef MATH_POINT_H
#define MATH_POINT_H

#include <Eigen/Core>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<int N, typename T>
class Direction;

template<int N, typename T>
class Vector;

template<int N, typename T>
class Point;

template<int N, typename T>
Vector<N, T> operator-(Point<N, T> const& to, Point<N, T> const& from);

template<int N, typename T>
bool operator!=(Point<N, T> const& p1, Point<N, T> const& p2);

template<int N, typename T = double>
class Point
{
 public:
  using scalar_type = T;
  using eigen_type = Eigen::Matrix<T, N, 1>;
  using vector_type = Vector<N, T>;
  using direction_type = Direction<N, T>;

 private:
  eigen_type p_;

 public:
  Point() = default;
  Point(eigen_type const& p) : p_(p) { }
  template<typename... U>
    requires (sizeof...(U) == N &&                                              // Exactly N coefficients.
      (std::convertible_to<U, T> && ...) &&                                     // All convertible to T.
      !(N == 1 && ((std::same_as<std::remove_cvref_t<U>, Point> || ...))))      // Do not replace copy/move constructor.
  constexpr Point(U&&... xs) : p_(static_cast<T>(std::forward<U>(xs))...) { }

  eigen_type& eigen() { return p_; }
  eigen_type const& eigen() const { return p_; }

  T& operator[](int i) { return p_(i); }
  T const& operator[](int i) const { return p_(i); }
  T& x() { static_assert(N >= 1); return p_[0]; }
  T& y() { static_assert(N >= 2); return p_[1]; }
  T& z() { static_assert(N >= 3); return p_[2]; }
  T x() const { static_assert(N >= 1); return p_[0]; }
  T y() const { static_assert(N >= 2); return p_[1]; }
  T z() const { static_assert(N >= 3); return p_[2]; }

  Point operator+(direction_type const& direction);
  Point operator+(vector_type const& v);
  Point operator-(vector_type const& v);
  Point& operator+=(direction_type const& direction);
  Point& operator+=(vector_type const& v);
  Point& operator-=(vector_type const& v);

  template<int N2, typename T2>
  friend Vector<N2, T2> operator-(Point<N2, T2> const& to, Point<N2, T2> const& from);

  template<int N2, typename T2>
  friend bool operator!=(Point<N2, T2> const& p1, Point<N2, T2> const& p2);

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

} // namespace math

#endif // MATH_POINT_H

#ifndef MATH_DIRECTION_H
#include "Direction.h"
#endif

#ifndef MATH_VECTOR_H
#include "Vector.h"
#endif

#ifndef MATH_POINT_H_definitions
#define MATH_POINT_H_definitions

namespace math {

template<int N, typename T>
Point<N, T> Point<N, T>::operator+(Direction<N, T> const& direction)
{
  return {p_ + direction.eigen()};
}

template<int N, typename T>
Point<N, T> Point<N, T>::operator+(Vector<N, T> const& v)
{
  return {p_ + v.eigen()};
}

template<int N, typename T>
Point<N, T> Point<N, T>::operator-(Vector<N, T> const& v)
{
  return {p_ - v.eigen()};
}

template<int N, typename T>
Point<N, T>& Point<N, T>::operator+=(Direction<N, T> const& direction)
{
  p_ += direction.eigen();
  return *this;
}

template<int N, typename T>
Point<N, T>& Point<N, T>::operator+=(Vector<N, T> const& v)
{
  p_ += v.eigen();
  return *this;
}

template<int N, typename T>
Point<N, T>& Point<N, T>::operator-=(Vector<N, T> const& v)
{
  p_ -= v.eigen();
  return *this;
}

template<int N, typename T>
Vector<N, T> operator-(Point<N, T> const& to, Point<N, T> const& from)
{
  return {from, to};
}

template<int N, typename T>
bool operator!=(Point<N, T> const& p1, Point<N, T> const& p2)
{
  for (int i = 0; i < N; ++i)
    if (p1.p_(i) != p2.p_(i))
      return true;
  return false;
}

#ifdef CWDEBUG
template<int N, typename T>
void Point<N, T>::print_on(std::ostream& os) const
{
  os << "{p_:" << p_ << '}';
}
#endif

} // namespace math

#endif // MATH_POINT_H_definitions
