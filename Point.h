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
  using eigen_type = Eigen::Matrix<T, N, 1>;

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

template<typename DerivedTypes>
struct PointOps
{
  static constexpr int derived_n = DerivedTypes::n;
  using derived_scalar_type    = typename DerivedTypes::scalar_type;
  using derived_type           = typename DerivedTypes::derived_type;

 private:
  auto& raw_() { return static_cast<derived_type*>(this)->raw(); }
  auto const& raw_() const { return static_cast<derived_type const*>(this)->raw(); }

 public:
  derived_scalar_type&       operator[](int i)       { return raw_().operator[](i); }
  derived_scalar_type const& operator[](int i) const { return raw_().operator[](i); }
  derived_scalar_type& x()      { return raw_().x(); }
  derived_scalar_type x() const { return raw_().x(); }
  derived_scalar_type& y()      { return raw_().y(); }
  derived_scalar_type y() const { return raw_().y(); }
  derived_scalar_type& z()      { return raw_().z(); }
  derived_scalar_type z() const { return raw_().z(); }

  inline derived_type operator+(typename DerivedTypes::direction_type const& direction) const;
  inline derived_type operator+(typename DerivedTypes::vector_type const& v) const;
  inline derived_type operator-(typename DerivedTypes::vector_type const& v) const;
  inline derived_type& operator+=(typename DerivedTypes::direction_type const& direction);
  inline derived_type& operator+=(typename DerivedTypes::vector_type const& v);
  inline derived_type& operator-=(typename DerivedTypes::vector_type const& v);
};

// Forward declaration, required for PointTypes<N, T>::derived_type.
template<int N, typename T>
class Point;

template<int N, typename T>
struct PointTypes
{
  static constexpr int n = N;
  using scalar_type    = T;
  // The following types will be incomplete while begin used in PointOps. Used for casting and/or function prototypes.
  using derived_type   = Point<N, T>;
  using vector_type    = Vector<N, T>;
  using direction_type = Direction<N, T>;
};

// Specialization of the Point operators specifically for math::Point itself.
template<int N, typename T>
struct PointOps<PointTypes<N, T>>
{
 private:
  // Get underlying Eigen type.
  auto& eigen_() { return static_cast<Point<N, T>*>(this)->eigen(); }
  auto const& eigen_() const { return static_cast<Point<N, T> const*>(this)->eigen(); }

 public:
  // Access coordinates.
  T&       operator[](int i)       { return eigen_()(i); }
  T const& operator[](int i) const { return eigen_()(i); }
  T& x()      { static_assert(N >= 1); return eigen_()[0]; }
  T x() const { static_assert(N >= 1); return eigen_()[0]; }
  T& y()      { static_assert(N >= 2); return eigen_()[1]; }
  T y() const { static_assert(N >= 2); return eigen_()[1]; }
  T& z()      { static_assert(N >= 3); return eigen_()[2]; }
  T z() const { static_assert(N >= 3); return eigen_()[2]; }

  // Add a direction.
  inline Point<N, T> operator+(Direction<N, T> const& direction) const;
  inline Point<N, T>& operator+=(Direction<N, T> const& direction);

  // Add a vector.
  inline Point<N, T> operator+(Vector<N, T> const& v) const;
  inline Point<N, T>& operator+=(Vector<N, T> const& v);

  // Subtract a vector.
  inline Point<N, T> operator-(Vector<N, T> const& v) const;
  inline Point<N, T>& operator-=(Vector<N, T> const& v);
};

template<int N, typename T = double>
class Point :
  public PointData<N, T>,               // Wraps underlying Eigen type and defines constructors.
  public PointOps<PointTypes<N, T>>     // Defines possible operations on the Point (uses the above specialization).
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

template<typename DerivedTypes>
typename DerivedTypes::derived_type PointOps<DerivedTypes>::operator+(typename DerivedTypes::direction_type const& direction) const
{
  return static_cast<derived_type>(raw_().operator+(direction.raw()));
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type PointOps<DerivedTypes>::operator+(typename DerivedTypes::vector_type const& v) const
{
  return static_cast<derived_type>(raw_().operator+(v.raw()));
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type PointOps<DerivedTypes>::operator-(typename DerivedTypes::vector_type const& v) const
{
  return static_cast<derived_type>(raw_().operator-(v.raw()));
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type& PointOps<DerivedTypes>::operator+=(typename DerivedTypes::direction_type const& direction)
{
  return static_cast<derived_type&>(raw_().operator+=(direction.raw()));
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type& PointOps<DerivedTypes>::operator+=(typename DerivedTypes::vector_type const& v)
{
  return static_cast<derived_type&>(raw_().operator+=(v.raw()));
}

template<typename DerivedTypes>
typename DerivedTypes::derived_type& PointOps<DerivedTypes>::operator-=(typename DerivedTypes::vector_type const& v)
{
  return static_cast<derived_type&>(raw_().operator-=(v.raw()));
}

// Specialization of the Point binary operators specifically for math::Point itself.

template<int N, typename T>
Point<N, T> PointOps<PointTypes<N, T>>::operator+(Direction<N, T> const& direction) const
{
  return {eigen_() + direction.eigen()};
}

template<int N, typename T>
Point<N, T> PointOps<PointTypes<N, T>>::operator+(Vector<N, T> const& v) const
{
  return {eigen_() + v.eigen()};
}

template<int N, typename T>
Point<N, T> PointOps<PointTypes<N, T>>::operator-(Vector<N, T> const& v) const
{
  return {eigen_() - v.eigen()};
}

template<int N, typename T>
Point<N, T>& PointOps<PointTypes<N, T>>::operator+=(Direction<N, T> const& direction)
{
  eigen_() += direction.eigen();
  return static_cast<Point<N, T>&>(*this);
}

template<int N, typename T>
Point<N, T>& PointOps<PointTypes<N, T>>::operator+=(Vector<N, T> const& v)
{
  eigen_() += v.eigen();
  return static_cast<Point<N, T>&>(*this);
}

template<int N, typename T>
Point<N, T>& PointOps<PointTypes<N, T>>::operator-=(Vector<N, T> const& v)
{
  eigen_() -= v.eigen();
  return static_cast<Point<N, T>&>(*this);
}

//-----------------------------------------------------------------------------
// Free function binary operators.
//

template<typename DerivedTypes>
typename DerivedTypes::vector_type operator-(PointOps<DerivedTypes> const& to, PointOps<DerivedTypes> const& from)
{
  return {static_cast<typename DerivedTypes::derived_type const&>(from), static_cast<typename DerivedTypes::derived_type const&>(to)};
}

template<typename DerivedTypes>
bool operator!=(PointOps<DerivedTypes> const& p1, PointOps<DerivedTypes> const& p2)
{
  for (int i = 0; i < DerivedTypes::n; ++i)
    if (p1[i] != p2[i])
      return true;
  return false;
}

//
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
