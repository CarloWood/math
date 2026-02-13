//NOT RECOVERED:
#pragma once

#include "TranslationVector.h"
#include "Direction.h"
#include "Vector.h"
#include <cmath>
#include <concepts>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <utility>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
// Transform defines a print_on method.
using utils::has_print_on::operator<<;
#endif

//------------------------------------------------------------------------------
// Affine transforms
//
// Transform used to be implemented directly with Qt's QTransform. In order to
// keep cairowindow usable without Qt, Transform is now parameterized with an
// affine backend that supplies the required operations.
//
// The backend operates on doubles. Qt projects can provide an adapter type that
// wraps QTransform and satisfies AffineTransformConcept.
//
template<typename AffineTransform>
concept AffineTransformConcept =
  std::default_initializable<AffineTransform> &&
  std::copy_constructible<AffineTransform> &&
  requires(AffineTransform transform, AffineTransform const const_transform,
           double x, double y, double dx, double dy, double factor, double radians)
  {
    { transform.m11() } -> std::same_as<double>;
    { transform.m21() } -> std::same_as<double>;
    { transform.m31() } -> std::same_as<double>;
    { transform.m12() } -> std::same_as<double>;
    { transform.m22() } -> std::same_as<double>;
    { transform.m32() } -> std::same_as<double>;
    { transform.translate(dx, dy) } -> std::same_as<AffineTransform&>;
    { transform.scale(factor, factor) } -> std::same_as<AffineTransform&>;
    { transform.rotate(radians) } -> std::same_as<AffineTransform&>;
    { const_transform.inverted() } -> std::same_as<AffineTransform>;
    { const_transform.map_point(x, y) } -> std::same_as<std::pair<double, double>>;
    { const_transform.map_vector(dx, dy) } -> std::same_as<std::pair<double, double>>;
    { const_transform * const_transform } -> std::same_as<AffineTransform>;
  };

// Standalone Affine backend.
//
// Matrix layout matches QTransform:
//
//   ⎡ m₁₁ m₁₂ 0 ⎤
//   ⎢ m₂₁ m₂₂ 0 ⎥
//   ⎣ m₃₁ m₃₂ 1 ⎦
//
// Mapping is done as a row-vector multiplication from the left:
//   [x y 1] * M
//
class AffineTransform2D
{
 private:
  // Identity.
  double m11_ = 1.0, m12_ = 0.0;
  double m21_ = 0.0, m22_ = 1.0;
  double m31_ = 0.0, m32_ = 0.0;

 private:
  static AffineTransform2D translation(double dx, double dy)
  {
    AffineTransform2D t;
    t.m31_ = dx;
    t.m32_ = dy;
    return t;
  }

  static AffineTransform2D scaling(double x_factor, double y_factor)
  {
    AffineTransform2D s;
    s.m11_ = x_factor;
    s.m22_ = y_factor;
    return s;
  }

  static AffineTransform2D rotation(double radians)
  {
    AffineTransform2D r;
    double const c = std::cos(radians);
    double const s = std::sin(radians);
    r.m11_ = c;
    r.m12_ = s;
    r.m21_ = -s;
    r.m22_ = c;
    return r;
  }

 public:
  AffineTransform2D() = default;

  // Matrix element accessors.
  double m11() const { return m11_; }
  double m21() const { return m21_; }
  double m31() const { return m31_; }
  double m12() const { return m12_; }
  double m22() const { return m22_; }
  double m32() const { return m32_; }

  // Apply a translation *before* the current transform.
  //
  //   returns
  //
  //     ⎡  1  0  0 ⎤⎡ m₁₁ m₁₂ 0 ⎤   ⎡ m₁₁                  m₁₂                  0 ⎤
  //     ⎢  0  1  0 ⎥⎢ m₂₁ m₂₂ 0 ⎥ = ⎢ m₂₁                  m₂₂                  0 ⎥
  //     ⎣ dx dy  1 ⎦⎣ m₃₁ m₃₂ 1 ⎦   ⎣ m₁₁dx + m₂₁dy + m₃₁  m₁₂dx + m₂₂dy + m₃₂  1 ⎦
  //
  AffineTransform2D& translate(double dx, double dy)
  {
    *this = translation(dx, dy) * *this;
    return *this;
  }

  // Apply a scaling *before* the current transform.
  //
  // - scale(s)
  //
  //   returns
  //
  //     ⎡ sx  0   0 ⎤⎡ m₁₁ m₁₂ 0 ⎤   ⎡ sx * m₁₁  sx * m₁₂  0 ⎤
  //     ⎢ 0   sy  0 ⎥⎢ m₂₁ m₂₂ 0 ⎥ = ⎢ sy * m₂₁  sy * m₂₂  0 ⎥
  //     ⎣ 0   0   1 ⎦⎣ m₃₁ m₃₂ 1 ⎦   ⎣      m₃₁       m₃₂  1 ⎦
  //
  AffineTransform2D& scale(double x_factor, double y_factor)
  {
    *this = scaling(x_factor, y_factor) * *this;
    return *this;
  }

  // Apply a rotation *before* the current transform.
  //
  // - rotate(radians)
  //
  //   returns
  //
  //     ⎡ cos(θ) sin(θ)  0 ⎤⎡ m₁₁ m₁₂ 0 ⎤   ⎡  cos(θ) * m₁₁ + sin(θ) * m₂₁   cos(θ) * m₁₂ + sin(θ) * m₂₂  0 ⎤
  //     ⎢-sin(θ) cos(θ)  0 ⎥⎢ m₂₁ m₂₂ 0 ⎥ = ⎢ -sin(θ) * m₁₁ + cos(θ) * m₂₁  -sin(θ) * m₁₂ + cos(θ) * m₂₂  0 ⎥
  //     ⎣ 0      0       1 ⎦⎣ m₃₁ m₃₂ 1 ⎦   ⎣  m₃₁                           m₃₂                          1 ⎦
  //
  AffineTransform2D& rotate(double radians)
  {
    *this = rotation(radians) * *this;
    return *this;
  }

  // - map_point(x, y)
  //
  //   returns
  //
  //           ⎡ m₁₁ m₁₂ 0 ⎤
  //    [x y 1]⎢ m₂₁ m₂₂ 0 ⎥ = [ x * m₁₁ + y * m₂₁ + m₃₁  x * m₁₂ + y * m₂₂ + m₃₂  1 ]
  //           ⎣ m₃₁ m₃₂ 1 ⎦
  //
  [[nodiscard]] std::pair<double, double> map_point(double x, double y) const
  {
    return {x * m11_ + y * m21_ + m31_, x * m12_ + y * m22_ + m32_};
  }

  [[nodiscard]] Point<2> map_point(Point<2> const& p) const
  {
    return {p.x() * m11_ + p.y() * m21_ + m31_, p.x() * m12_ + p.y() * m22_ + m32_};
  }

  // Same, without translation (assume m₃₁ = m₃₂ = 0).
  [[nodiscard]] std::pair<double, double> map_vector(double dx, double dy) const
  {
    return {dx * m11_ + dy * m21_, dx * m12_ + dy * m22_};
  }

  [[nodiscard]] Vector<2> map_vector(Vector<2> const& v) const
  {
    return {v.x() * m11_ + v.y() * m21_, v.x() * m12_ + v.y() * m22_};
  }

  [[nodiscard]] Direction<2> map_direction(Direction<2> const& v) const
  {
    return Direction<2>{Direction<2>::eigen_type{v.x() * m11_ + v.y() * m21_, v.x() * m12_ + v.y() * m22_}};
  }

  // - inverted()
  //
  //   returns
  //
  //            ⎡  m₂₂           -m₁₂            0   ⎤
  //    (1/det) ⎢ -m₂₁            m₁₁            0   ⎥
  //            ⎣  m₂₁m₃₂-m₂₂m₃₁  m₁₂m₃₁-m₁₁m₃₂  det ⎦
  //
  //   where det = m₁₁m₂₂-m₁₂m₂₁
  //
  [[nodiscard]] AffineTransform2D inverted() const
  {
    double const det = m11_ * m22_ - m12_ * m21_;
    // Don't scale with a factor of zero.
    ASSERT(det != 0.0);
    double const inv_det = 1.0 / det;

    AffineTransform2D inv;
    inv.m11_ =  m22_ * inv_det;
    inv.m12_ = -m12_ * inv_det;
    inv.m21_ = -m21_ * inv_det;
    inv.m22_ =  m11_ * inv_det;
    inv.m31_ = -(inv.m21_ * m32_ + inv.m11_ * m31_);    // (m₂₁m₃₂-m₂₂m₃₁)/det = -((-m₂₁/det) m₃₂ + (m₂₂/det) m₃₁)
    inv.m32_ = -(inv.m12_ * m31_ + inv.m22_ * m32_);    // (m₁₂m₃₁-m₁₁m₃₂)/det = -((-m₁₂/det) m₃₁ + (m₁₁/det) m₃₂)
    return inv;
  }

  //
  //     ⎡ l₁₁ l₁₂ 0 ⎤⎡ r₁₁ r₁₂ 0 ⎤   ⎡ l₁₁r₁₁+l₁₂r₂₁      l₁₁r₁₂+l₁₂r₂₂      0 ⎤
  //     ⎢ l₂₁ l₂₂ 0 ⎥⎢ r₂₁ r₂₂ 0 ⎥ = ⎢ l₂₁r₁₁+l₂₂r₂₁      l₂₁r₁₂+l₂₂r₂₂      0 ⎥
  //     ⎣ l₃₁ l₃₂ 1 ⎦⎣ r₃₁ r₃₂ 1 ⎦   ⎣ l₃₁r₁₁+l₃₂r₂₁+r₃₁  l₃₁r₁₂+l₃₂r₂₂+r₃₂  1 ⎦
  //
  friend AffineTransform2D operator*(AffineTransform2D const& lhs, AffineTransform2D const& rhs)
  {
    AffineTransform2D result;
    result.m11_ = lhs.m11_ * rhs.m11_ + lhs.m12_ * rhs.m21_;
    result.m12_ = lhs.m11_ * rhs.m12_ + lhs.m12_ * rhs.m22_;
    result.m21_ = lhs.m21_ * rhs.m11_ + lhs.m22_ * rhs.m21_;
    result.m22_ = lhs.m21_ * rhs.m12_ + lhs.m22_ * rhs.m22_;
    result.m31_ = lhs.m31_ * rhs.m11_ + lhs.m32_ * rhs.m21_ + rhs.m31_;
    result.m32_ = lhs.m31_ * rhs.m12_ + lhs.m32_ * rhs.m22_ + rhs.m32_;
// RECOVERED:
    return result;
  }
};

// Verify that AffineTransform2D fullfills the concept.
static_assert(AffineTransformConcept<AffineTransform2D>);

template<CS from_cs, CS to_cs, bool inverted_ = false, AffineTransformConcept AffineTransformBackend = AffineTransform2D>
class Transform
{
 private:
  template<CS from_cs2, CS to_cs2, bool inverted_2, AffineTransformConcept AffineTransformBackend2>
  friend class Transform;

  AffineTransformBackend m_;

 private:
  Transform(AffineTransformBackend const& m) : m_(m) { }

 public:
  // Construct an Identity Transform.
  Transform() = default;

  // Prepend an affine transform to the current Transform.
  Transform& translate(TranslationVector<to_cs> const& tv) requires (inverted_ == false);
  Transform& scale(double x_factor, double y_factor) requires (inverted_ == false);
  Transform& scale(double factor) { return scale(factor, factor); }
  Transform& rotate(double radians) requires (inverted_ == false);

  // Map point (x,y).
  [[nodiscard]] std::pair<double, double> map_point(double x, double y) const
  requires (inverted_ == false)
  {
    return m_.map_point(x, y);
  }

  [[nodiscard]] cs::Point<to_cs> map_point(cs::Point<from_cs> const& p) const
  requires (inverted_ == false)
  {
    return cs::Point<to_cs>{m_.map_point(p.raw())};
  }

  [[nodiscard]] std::pair<double, double> map_vector(double dx, double dy) const
  requires (inverted_ == false)
  {
    return m_.map_vector(dx, dy);
  }

  [[nodiscard]] cs::Vector<to_cs> map_vector(cs::Vector<from_cs> const& v) const
  requires (inverted_ == false)
  {
    return cs::Vector<to_cs>{m_.map_vector(v.raw())};
  }

  [[nodiscard]] cs::Direction<to_cs> map_direction(cs::Direction<from_cs> const& d) const
  {
    return cs::Direction<to_cs>{m_.map_direction(d.raw())};
  }

  // Just scale.
  //
  [[nodiscard]] cs::Size<to_cs> map_size(cs::Size<from_cs> const& s) const
  {
    double const x_factor = std::hypot(m_.m11(), m_.m12());
    double const y_factor = std::hypot(m_.m21(), m_.m22());

    // Just scale; scale the X and Y axis vectors by the full linear part.
    if constexpr (!inverted_)
      return {s.width() * x_factor, s.height() * y_factor};
    else
      return {s.width() / x_factor, s.height() / y_factor};
  }

  // The inverted Transfrom converts from `to_cs` to `from_cs`.
  Transform<to_cs, from_cs, false, AffineTransformBackend> inverted() const
  {
    return {inverted_ ? m_ : m_.inverted()};
  }

  TranslationVector<to_cs> translation() const
  requires (inverted_ == false)
  {
    return {m_.m31(), m_.m32()};
  }

  double x_scale() const
  requires (inverted_ == false)
  {
    double const x_factor = std::hypot(m_.m11(), m_.m12());
    return x_factor;
  }

// NOT RECOVERED:
  double y_scale() const
  {
    double const y_factor = std::hypot(m_.m21(), m_.m22());
    return y_factor;
  }

  std::pair<double, double> scale_factors() const
  {
    return {x_scale(), y_scale()};
  };


  // Return the direction of the mapped x-axis.
  //
  // This is the direction of the basis vector (1, 0) after applying the linear
  // part of this transform; it is well-defined even when the transform contains
  // non-uniform scaling and/or shear.
  cs::Direction<to_cs> x_axis_direction() const
  {
    return cs::Direction<to_cs>(Point<2>(m_.m11(), m_.m12()));
  }

  // Return the direction of the mapped y-axis (basis vector (0, 1) mapped by the linear part).
  cs::Direction<to_cs> y_axis_direction() const
  {
    return cs::Direction<to_cs>(Point<2>(m_.m21(), m_.m22()));
  }

  // The inverse converts from `to_cs` to `from_cs`!
  Transform<to_cs, from_cs, !inverted_, AffineTransformBackend> const& inverse() const
  {
    return reinterpret_cast<Transform<to_cs, from_cs, !inverted_, AffineTransformBackend> const&>(*this);
  }

  // Let A_M1_B be non-inverted and convert from A to B.
  // Let B_M2_C be non-inverted and convert from B to C.
  //
  // Multiplication between two non-inverted Transforms.
  // 1. A_M12_C = A_M1_B * B_M2_C
  //
  //
  // Let A_M1_B^-1 be an inverted matrix that converts from B to A, and therefore denote it as B_M1inv_A
  // Let B_M2_C^-1 be an inverted matrix that converts from C to B, and therefore denote it as C_M2inv_B
  //
  // Multiplication between two inverted Transforms.
  // 2. A_M34inv_C = A_M4inv_B * B_M3inv_C
  //
  // Note that C_M34_A = C_M3_B * B_M4_A
  //
  //
  // Let A_M5inv_B = B_M5_A^-1 be an inverted matrix that converts from A to B.
  // Let B_M56_C = B_M5_A * A_M6_C.
  //
  // Multiplication between an inverted Transform and a non-inverted Transform.
  // 3. A_M6_C = A_M5inv_B * B_M56_C
  //
  //
  // Let A_M78_B = A_M7_C * C_M8_B.
  //
  // Multiplication between a non-inverted Transform and an inverted Transform.
  // 4. A_M7_C = A_M78_B * B_M8inv_C
  //
  //
  // Then using specializations, where from_cs = A, to_cs = B and result_cs = C we'd have:
  //
  // Specialization for 1 (neither input inverted, result also not inverted).
  // std::enable_if_t<!inverted_, Transform<from_cs, result_cs, false>> operator*(Transform<to_cs, result_cs, false> const& rhs) const;
  //
  // Specialization for 2 (both inputs inverted, output also inverted).
  // std::enable_if_t<inverted_, Transform<from_cs, result_cs, true>> operator*(Transform<to_cs, result_cs, true> const& rhs) const;
  //
  // Specialization for 3 (only lhs input inverted, output not inverted).
  // std::enable_if_t<inverted_, Transform<from_cs, result_cs, false>> operator*(Transform<to_cs, result_cs, false> const& rhs) const;
  //
  // Specialization for 4 (only rhs input inverted, output not inverted).
  //std::enable_if_t<!inverted_, Transform<from_cs, result_cs, false>> operator*(Transform<to_cs, result_cs, true> const& rhs) const;
  //
  template<CS result_cs, bool rhs_inverted>
  Transform<from_cs, result_cs, inverted_ && rhs_inverted, AffineTransformBackend>
  operator*(Transform<to_cs, result_cs, rhs_inverted, AffineTransformBackend> const& rhs) const
  {
    // 1. Multiplication between two non-inverted Transforms.
    if constexpr (!inverted_ && !rhs_inverted)
    {
      return {m_ * rhs.m_};
    }
    // 2. Multiplication between two inverted Transforms.
    else if constexpr (inverted_ && rhs_inverted)
    {
      // A^-1 * B^-1 = (B * A)^-1
      return {rhs.m_ * m_};
    }
    // 3. Multiplication between an inverted Transform and a non-inverted Transform.
    else if constexpr (inverted_ && !rhs_inverted)
    {
      return {m_.inverted() * rhs.m_};
    }
    // 4. Multiplication between a non-inverted Transform and an inverted Transform.
    else if constexpr (!inverted_ && rhs_inverted)
    {
      return {m_ * rhs.m_.inverted()};
    }
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    std::ostringstream prefix;
    prefix << utils::to_string(from_cs) << "_transform_" << utils::to_string(to_cs) << ":";
    int const prefix_len = std::max((int)prefix.str().length(), 24);

    os << '\n' << std::setw(prefix_len) << " " <<
      std::left << "⎛" << std::setw(9) << m_.m11() << " " << std::setw(9) << m_.m12() << " " << 0.0 << "⎞" << std::right;
    os << '\n' << std::right << std::setw(prefix_len) << prefix.str() <<
      std::left << "⎜" << std::setw(9) << m_.m21() << " " << std::setw(9) << m_.m22() << " " << 0.0 << "⎟" << std::right;
    os << '\n' << std::setw(prefix_len) << " " <<
      std::left << "⎝" << std::setw(9) << m_.m31()  << " " << std::setw(9) << m_.m32() << " " << 1.0 << "⎠" << std::right;
  }
#endif
};

template<CS from_cs, CS to_cs, bool inverted_, AffineTransformConcept AffineTransformBackend>
Transform<from_cs, to_cs, inverted_, AffineTransformBackend>&
Transform<from_cs, to_cs, inverted_, AffineTransformBackend>::translate(TranslationVector<to_cs> const& tv)
requires (inverted_ == false)
{
  m_.translate(tv.as_vector().x(), tv.as_vector().y());
  return *this;
}

template<CS from_cs, CS to_cs, bool inverted_, AffineTransformConcept AffineTransformBackend>
Transform<from_cs, to_cs, inverted_, AffineTransformBackend>&
Transform<from_cs, to_cs, inverted_, AffineTransformBackend>::scale(double x_factor, double y_factor)
requires (inverted_ == false)
{
  m_.scale(x_factor, y_factor);
  return *this;
}

template<CS from_cs, CS to_cs, bool inverted_, AffineTransformConcept AffineTransformBackend>
Transform<from_cs, to_cs, inverted_, AffineTransformBackend>&
Transform<from_cs, to_cs, inverted_, AffineTransformBackend>::rotate(double radians)
requires (inverted_ == false)
{
  m_.rotate(radians);
  return *this;
}

} // namespace math
