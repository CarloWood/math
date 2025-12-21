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
class Point;

template<int N, typename T>
class Vector;

template<int N, typename T>
class Direction;

template<int N, typename T>
class PointData
{
 public:
  static constexpr int n = N;
  using scalar_type = T;
  using eigen_type = Eigen::Matrix<T, N, 1>;
  using vector_type = Vector<N, T>;
  using direction_type = Direction<N, T>;

 protected:
  eigen_type p_;

 public:
  PointData() = default;
  PointData(eigen_type const& p) : p_(p) { }
  template<typename... U>
    requires (sizeof...(U) == N &&                                              // Exactly N coefficients.
      (std::convertible_to<U, T> && ...) &&                                     // All convertible to T.
      !(N == 1 && ((std::same_as<std::remove_cvref_t<U>, PointData> || ...))))  // Do not replace copy/move constructor.
  constexpr PointData(U&&... xs) : p_(static_cast<T>(std::forward<U>(xs))...) { }

  eigen_type& eigen() { return p_; }
  eigen_type const& eigen() const { return p_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

template<typename Derived, typename Base>
struct PointOps
{
  static constexpr int derived_n = Base::n;
  using derived_scalar_type = typename Base::scalar_type;
  using vector_type = typename Base::vector_type;
  using direction_type = typename Base::direction_type;
  using raw_type = Point<Base::n, derived_scalar_type>;

  static constexpr bool is_raw = std::is_same_v<Derived, raw_type>;

  derived_scalar_type& operator[](int i)
  {
    if constexpr (is_raw)
      return p()(i);
    else
      return r().operator[](i);
  }

  derived_scalar_type const& operator[](int i) const
  {
    if constexpr (is_raw)
      return p()(i);
    else
      return r().operator[](i);
  }

  derived_scalar_type& x()
  {
    if constexpr (is_raw)
    {
      static_assert(derived_n >= 1);
      return p()[0];
    }
    else
      return r().x();
  }

  derived_scalar_type& y()
  {
    if constexpr (is_raw)
    {
      static_assert(derived_n >= 2);
      return p()[1];
    }
    else
      return r().y();
  }

  derived_scalar_type& z()
  {
    if constexpr (is_raw)
    {
      static_assert(derived_n >= 3);
      return p()[2];
    }
    else
      return r().z();
  }

  derived_scalar_type x() const
  {
    if constexpr (is_raw)
    {
      static_assert(derived_n >= 1);
      return p()[0];
    }
    else
      return r().x();
  }

  derived_scalar_type y() const
  {
    if constexpr (is_raw)
    {
      static_assert(derived_n >= 2);
      return p()[1];
    }
    else
      return r().y();
  }

  derived_scalar_type z() const
  {
    if constexpr (is_raw)
    {
      static_assert(derived_n >= 3);
      return p()[2];
    }
    else
      return r().z();
  }

  Derived operator+(direction_type const& direction);
  Derived operator+(vector_type const& v);
  Derived operator-(vector_type const& v);
  Derived& operator+=(direction_type const& direction);
  Derived& operator+=(vector_type const& v);
  Derived& operator-=(vector_type const& v);

 private:
  auto& p() { return static_cast<Derived*>(this)->eigen(); }
  auto const& p() const { return static_cast<Derived const*>(this)->eigen(); }
  auto& r() { return static_cast<Derived*>(this)->raw(); }
  auto const& r() const { return static_cast<Derived const*>(this)->raw(); }
};

template<int N, typename T = double>
class Point : public PointData<N, T>, public PointOps<Point<N, T>, PointData<N, T>>
{
 public:
  // Constructors.
  using PointData<N, T>::PointData;
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

template<typename Derived, typename Base>
Derived PointOps<Derived, Base>::operator+(direction_type const& direction)
{
  if constexpr (is_raw)
    return static_cast<Derived>(p() + direction.eigen());
  else
    return static_cast<Derived>(r().operator+(direction));
}

template<typename Derived, typename Base>
Derived PointOps<Derived, Base>::operator+(vector_type const& v)
{
  if constexpr (is_raw)
    return static_cast<Derived>(p() + v.eigen());
  else
    return static_cast<Derived>(r().operator+(v));
}

template<typename Derived, typename Base>
Derived PointOps<Derived, Base>::operator-(vector_type const& v)
{
  if constexpr (is_raw)
    return static_cast<Derived>(p() - v.eigen());
  else
    return static_cast<Derived>(r().operator-(v));
}

template<typename Derived, typename Base>
Derived& PointOps<Derived, Base>::operator+=(direction_type const& direction)
{
  if constexpr (is_raw)
    return static_cast<Derived&>(p() += direction.eigen());
  else
    return static_cast<Derived&>(r().operator+=(direction));
}

template<typename Derived, typename Base>
Derived& PointOps<Derived, Base>::operator+=(vector_type const& v)
{
  if constexpr (is_raw)
    return static_cast<Derived&>(p() += v.eigen());
  else
    return static_cast<Derived&>(r().operator+=(v));
}

template<typename Derived, typename Base>
Derived& PointOps<Derived, Base>::operator-=(vector_type const& v)
{
  if constexpr (is_raw)
    return static_cast<Derived&>(p() -= v.eigen());
  else
    return static_cast<Derived&>(r().operator-=(v));
}

//-----------------------------------------------------------------------------
// Free function binary operators.
//
template<typename Derived, typename Base>
requires std::derived_from<Derived, PointOps<Derived, Base>>
typename Derived::vector_type operator-(PointOps<Derived, Base> const& to, PointOps<Derived, Base> const& from)
{
  return {static_cast<Derived const&>(from), static_cast<Derived const&>(to)};
}

template<typename Derived, typename Base>
requires std::derived_from<Derived, PointOps<Derived, Base>>
bool operator!=(PointOps<Derived, Base> const& p1, PointOps<Derived, Base> const& p2)
{
  for (int i = 0; i < Derived::n; ++i)
    if (p1[i] != p2[i])
      return true;
  return false;
}
//-----------------------------------------------------------------------------

#ifdef CWDEBUG
template<int N, typename T>
void PointData<N, T>::print_on(std::ostream& os) const
{
  os << "{p_:(";
  char const* separator = "";
  for (int i = 0; i < N; ++i)
  {
    os << separator << p_(i);
    separator = ", ";
  }
  os << ")}";
}
#endif

} // namespace math

#endif // MATH_POINT_H_definitions
